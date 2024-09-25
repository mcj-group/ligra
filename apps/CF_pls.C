//This is a Ligra implementation of a gradient descent-based parallel
//algorithm for collaborative filtering (matrix factorization)
//described in the paper "GraphMat: High performance graph analytics
//made productive", VLDB 2015
//(https://github.com/narayanan2004/GraphMat/blob/master/src/SGD.cpp). The
//Ligra implementation is written by Yunming Zhang (and slightly
//modified by Julian Shun).

//The input to the program is a weighted bipartite graph between users
//and items, where the weights represent the rating a user gives to an
//item. The optional arguments to the program are as follows: "-K"
//specifies the dimension of the latent vector (default is 20),
//"-numiter" is the number of iterations of gradient descent to run
//(default is 5), "-step" is the step size in the algorithm (default
//is 0.00000035), "-lambda" is the regularization parameter (default
//is 0.001), and "-randInit" specifies that the latent vector should
//be initialized randomly (by default every entry is initialized to
//0.5).

#define WEIGHTED 1
#include "ligra.h"
#include "ligra_pls.h"
#include "swarm/api.h"

//uncomment the following line to compute the sum of squared errors per iteration
//#define COMPUTE_ERROR 1
//uncomment the following line to print out the sum of values in latent vector
//#define DEBUG 1


#ifdef COMPUTE_ERROR
double* squaredErrors;
#endif

#define TS_PER_ITERATION (2)

template <class vertex>
struct CF_Edge_F {
  double* latent_curr, *error;
  vertex* V;
  int K;
  CF_Edge_F(vertex* _V, double* _latent_curr, double* _error, int _K) :
    latent_curr(_latent_curr), error(_error), V(_V), K(_K) {}
  //updates latent vector based on neighbors' data
  inline bool update(uintE s, uintE d, intE edgeLen) const {
    double estimate = 0;
    long current_offset = K*d, ngh_offset = K*s;
    double *cur_latent =  &latent_curr[current_offset], *ngh_latent = &latent_curr[ngh_offset];
    for(int i = 0; i < K; i++){
      estimate +=  cur_latent[i]*ngh_latent[i];
    }
    double err = edgeLen - estimate;

#ifdef COMPUTE_ERROR
    squaredErrors[d] += err*err;
#endif

    double* cur_error = &error[current_offset];
    for (int i = 0; i < K; i++){
      cur_error[i] += ngh_latent[i]*err;
    }
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen) const {
    //not needed as we will always do pull based
    return update(s,d,edgeLen);
  }

  inline bool cond(intT d) const { return cond_true(d); }

  swarm::Hint hintCG(uintE d) const {
      return {swarm::Hint::cacheLine(&error[K*d]), EnqFlags::MAYSPEC};
  }

  swarm::Hint hintFG(uintE s, uintE d) const { return hintCG(d); }
};

struct CF_Vertex_F {
  double step, lambda;
  double *latent_curr, *error;
  int K;
  CF_Vertex_F(double _step, double _lambda, double* _latent_curr, double* _error, int _K) :
    step(_step), lambda(_lambda), latent_curr(_latent_curr), error(_error), K(_K) {}
  inline bool operator () (uintE i) const {
    for (int j = 0; j < K; j++){
      latent_curr[K*i + j] += step*(-lambda*latent_curr[K*i + j] + error[K*i + j]);
      error[K*i+j] = 0.0;
    }
#ifdef COMPUTE_ERROR
    squaredErrors[i] = 0;
#endif
    return 1;
  }
  swarm::Hint hint(uintE v) const {
    return {swarm::Hint::cacheLine(&error[K*v]), EnqFlags::MAYSPEC};
  }
};

// [mcj] OMG we need a timestamped sequential loop in the runtime
template <class vertex>
struct Iteration {
    graph<vertex>& GA;
    vertexSubset& Frontier;
    swarm::Timestamp lastTS;
    CF_Edge_F<vertex> edge_f;
    CF_Vertex_F vertex_f;

    void operator() (swarm::Timestamp ts) {
        if (ts >= lastTS) return;
        swarm::info("iteration %ld", (ts - 1) / TS_PER_ITERATION);

        //edgemap to accumulate error for each node
        pls_edgeMapDense<pbbs::empty>(ts, GA, Frontier, &edge_f);

        swarm::enqueueLambda([this] (swarm::Timestamp ts) {
            //vertexmap to update the latent vectors
            pls_vertexMapDense(ts, this->Frontier, &this->vertex_f);
            swarm::enqueueLambda(this, ts + 1, EnqFlags(NOHINT | CANTSPEC));
        }, ts + 1,
        // Use CANTSPEC to treat ts + 1 as a barrier (this has better
        // performance than swarm::serialize, as it does not incur an exception,
        // sending the whole system into exception mode).
        EnqFlags(NOHINT | CANTSPEC));
    }
};


template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  int K = P.getOptionIntValue("-K",20); //dimensions of the latent vector
  int numIter = P.getOptionIntValue("-numiter",5); //number of iterations
  double step = P.getOptionDoubleValue("-step",0.00000035); //step size
  double lambda = P.getOptionDoubleValue("-lambda",0.001); //lambda
  bool randInit = P.getOption("-randInit"); //pass flag for random initialization of latent vector
  //initialize latent vectors and errors
  const intE n = GA.n;
  double* latent_curr = newA(double, K*n);
  double* error = newA(double, K*n);
#ifdef COMPUTE_ERROR
  squaredErrors = newA(double,n);
#endif
  if(randInit) {
    srand(0);
    long seed = rand();
    parallel_for(uintE i = 0; i < n; i++){
#ifdef COMPUTE_ERROR
      squaredErrors[n] = 0;
#endif
      for (int j = 0; j < K; j++){
        latent_curr[i*K+j] = ((double)(seed+hashInt((uintE)i*K+j))/(double)UINT_E_MAX);
        error[i*K+j] = 0.0;
      }
    }
  } else {
    parallel_for(uintE i = 0; i < n; i++){
#ifdef COMPUTE_ERROR
      squaredErrors[n] = 0;
#endif
      for (int j = 0; j < K; j++){
        latent_curr[i*K+j] = 0.5; //default initial value of 0.5
        error[i*K+j] = 0.0;
      }
    }
  }

  bool* frontier = newA(bool,n);
  vertexSubset Frontier(n,n,frontier);
  swarm::fill<EnqFlags(NOHINT | MAYSPEC)>(frontier, frontier + n, true, 0ul);

  swarm::Timestamp lastTS = 1ul + numIter * TS_PER_ITERATION;
  Iteration<vertex> iteration = {GA, Frontier, lastTS,
      CF_Edge_F<vertex>(GA.V,latent_curr,error,K),
      CF_Vertex_F(step,lambda,latent_curr,error,K)
    };
  swarm::enqueueLambda(&iteration, 1ul, EnqFlags(NOHINT | CANTSPEC));
  swarm::run();

  Frontier.del(); free(latent_curr); free(error);
#ifdef COMPUTE_ERROR
  free(squaredErrors);
#endif
}
