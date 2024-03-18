// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Triangle counting code (assumes a symmetric graph, so pass the "-s"
// flag). This is not optimized (no ordering heuristic is used)--for
// optimized code, see "Multicore Triangle Computations Without
// Tuning", ICDE 2015. Currently only works with uncompressed graphs,
// and not with compressed graphs.
#include "ligra.h"
#include "quickSort.h"

#include "ligra_pls.h"
#include "swarm/api.h"
#include "swarm/numeric.h"
#include "swarm/cps.h"
#include <functional>

//assumes sorted neighbor lists
template <class vertex>
long countCommon(const vertex& A, const vertex& B, uintE a, uintE b) {
  uintT i=0,j=0,nA = A.getOutDegree(), nB = B.getOutDegree();
  const uintE* nghA = A.getOutNeighbors();
  const uintE* nghB = B.getOutNeighbors();
  long ans=0;
  while (i < nA && j < nB && nghA[i] < a && nghB[j] < b) { //count "directed" triangles
    if (nghA[i]==nghB[j]) i++, j++, ans++;
    else if (nghA[i] < nghB[j]) i++;
    else j++;
  }
  return ans;
}

template <>
long countCommon<compressedSymmetricVertex>(
        const compressedSymmetricVertex&,
        const compressedSymmetricVertex&,
        uintE, uintE) {
    std::cerr << "Unimplemented" << std::endl;
    std::abort();
}

template <>
long countCommon<compressedAsymmetricVertex>(
        const compressedAsymmetricVertex&,
        const compressedAsymmetricVertex&,
        uintE, uintE) {
    std::cerr << "Unimplemented" << std::endl;
    std::abort();
}


template <class vertex>
struct countF { //for edgeMap
  vertex* V;
  long* counts;
  countF(vertex* _V, long* _counts) : V(_V), counts(_counts) {}
  inline bool update (uintE s, uintE d) {
    if(s > d) //only count "directed" triangles
      writeAdd(&counts[s], countCommon<vertex>(V[s],V[d],s,d));
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d) {
    if (s > d) //only count "directed" triangles
#ifdef NONATOMIC_TASKS
      writeAdd(&counts[s], countCommon<vertex>(V[s],V[d],s,d));
#else
      // No CAS because Hints guarantee non-speculative correctness
      counts[s] += countCommon<vertex>(V[s],V[d],s,d);
#endif
    return 1;
  }
  inline bool cond (uintE d) const { return cond_true(d); } //does nothing

  swarm::Hint hintCG(uintE d) const {
      // The CG task centers around vertex 'd', but update()s counts at all its
      // neighbors 's' via writeAdd, which uses compare-and-swap.
      // NOHINT: a hint of 'd' isn't particularly useful (unverified)
      // MAYSPEC: speculation is unnecessary because of the compare-and-swap
      return EnqFlags(NOHINT | MAYSPEC);
  }

  swarm::Hint hintFG(uintE s, uintE d) const {
      // s: the FG task updates a count for vertex s
      // MAYSPEC: we're using cacheLine hints, and compare-and-swap
      return {swarm::Hint::cacheLine(&counts[s]), EnqFlags::MAYSPEC};
  }
};

struct intLT { bool operator () (uintT a, uintT b) { return a < b; }; };

template <class vertex>
struct initF { //for vertexMap to initial counts and sort neighbors for merging
  vertex* V;
  initF(vertex* _V) : V(_V) {}
  inline bool operator () (uintE i) {
    quickSort(V[i].getOutNeighbors(),V[i].getOutDegree(),intLT());
    return 1;
  }

  swarm::Hint hint(uintE) const {
    return EnqFlags(NOHINT | MAYSPEC); // The quickSorts are independent
  }
};


template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  uintT n = GA.n;
  long* counts = newA(long,n);
  bool* frontier = newA(bool,n);
  vertexSubset Frontier(n,n,frontier);

  initF<vertex> init_f(GA.V);
  countF<vertex> count_f(GA.V, counts);

  //frontier contains all vertices
  swarm::fill<EnqFlags(NOHINT | MAYSPEC)>(frontier, frontier + n, true, 0ul);
  swarm::fill<EnqFlags(NOHINT | MAYSPEC)>(counts, counts + n, 0, 1ul);

  pls_vertexMapDense(1ul, Frontier, &init_f);

  pls_cbegin(2ul, EnqFlags(NOHINT | CANTSPEC),
             [&count_f, &GA, &Frontier, &counts]);
  pls_edgeMapDense<pbbs::empty>(ts, GA, Frontier, &count_f);

  pls_cbegin(3ul, EnqFlags(NOHINT | CANTSPEC), [&GA, &counts]);

  swarm::reduce(counts, counts + GA.n, 0ul, std::plus<uint64_t>(), ts,
      [] (swarm::Timestamp ts, uint64_t sum) {
        std::cout << "triangle count = " << sum << std::endl;
      });

  pls_cend();
  pls_cend();

  swarm::run();
  Frontier.del(); free(counts);
}
