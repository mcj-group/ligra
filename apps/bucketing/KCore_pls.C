#include "ligra.h"
#include <vector>
#include "index_map.h"
#include "bucket.h"
#include "edgeMapReduce.h"
#include "swarm/api.h"
#include "swarm/scheduler.h"
using namespace swarm;

#ifdef COMPETITIVE_SCHEDULE
Scheduler<EnqFlags::UPDATEABLE, true> *s;
#else
Scheduler<EnqFlags::MAYSPEC, true> *s;
#endif
#ifdef ACCUMULATE
uintE* degrees;
#endif

template<class vertex> struct Update;

//function because functors don't seem to work properly with enqueue macros
template<class vertex>
static void callUpdate(Timestamp ts, Update<vertex> *u)
{
    (*u)(ts, s->current_object());
}

#ifdef ACCUMULATE
//generate a new timestamp based on current timestamp and 
//remaining degree.  As we get closer to removing this vertex,
//our timestamps get more precise.
static Timestamp newTS(Timestamp curTS, uintE v) {
    assert(degrees[v] >= curTS);
    if (degrees[v] <= curTS + 2) return degrees[v];
    return (curTS >> 1) + (degrees[v] >> 1);
}

//Track degree outside the scheduler.
//If degree is decremented below timestamp, reenqueue with a new timestamp
template<class vertex>
static void decrementDegree(Timestamp ts, Update<vertex> *u, uintE v) {
    if (degrees[v] <= ts) return;
    Timestamp oldTS = s->extract_ts(v);
    assert(oldTS <= degrees[v]);
    degrees[v]--;
    if (oldTS > degrees[v]) 
        absolute_enqueue(s, callUpdate<vertex>, newTS(ts, v), v, u);
}
#else
//function to call the relative update for fine-grained implementations
template<class vertex> 
static void decrementDegree(Timestamp ts, Update<vertex> *u, uintE v)
{
    //swarm::info("extracting ts for vertex %lu", v);
#ifndef COMPETITIVE_SCHEDULE
    if (s->extract_ts(v) > ts) //early exit simulates COMPETITIVE semantics w/out the flag
    {
#endif
   //     swarm::info("got ts for vertex %lu", v);
        relative_enqueue(s, callUpdate<vertex>, -1, v, u);
#ifndef COMPETITIVE_SCHEDULE
    }
#endif
//    swarm::info("done decrement");
}
#endif

//functor for the neighborhood op bc the graph is templated so I can't just use a global
template <class vertex>
struct Update
{
    graph<vertex>& GA;
    
    void operator() (Timestamp ts, uintE v)
    {
#ifndef COMPETITIVE_SCHEDULE
        if (s->extract_ts(v) < ts) return; //early exit simulates COMPETITIVE semantics w/out the flag
#endif
#ifdef ACCUMULATE
        //If our timestamp is less than our degree, then we overshot
        //so reenqueue with a new timestamp
        if (ts < degrees[v]) 
            absolute_enqueue(s, callUpdate<vertex>, newTS(ts, v), v, this);
        else
#endif 
            enqueue_all_progressive<swarm::max_children>(
                swarm::u64it(0), swarm::u64it(GA.V[v].getOutDegree()),
                [=] (Timestamp ts, uint64_t i) {
                    uintE ngh = GA.V[v].getOutNeighbor(i);
#ifdef COARSE_GRAIN
#ifndef COMPETITIVE_SCHEDULE
                    decrementDegree<vertex>(ts, this, ngh); },
#else
                    relative_enqueue(s, callUpdate<vertex>, -1, ngh, this); },
#endif
#else
                    enqueue(decrementDegree<vertex>, ts, ngh, this, ngh); },
#endif
                [=] (uint64_t) { return ts; },
                [] (uint64_t) { return EnqFlags::NOHINT; }); //no point in a hint unless graph is well-ordered/low degree 
    }
};

template <class vertex>
static void init(Timestamp, Update<vertex> *u, uintE v) {
    s->set_base_ts(v, u->GA.V[v].getOutDegree());
}

template <class vertex>
array_imap<uintE> KCore(graph<vertex>& GA, size_t num_buckets=128) {
  const size_t n = GA.n; const size_t m = GA.m;
  auto D = array_imap<uintE>(n, [&] (size_t i) { return GA.V[i].getOutDegree(); });
#ifdef ACCUMULATE
  degrees = D.s;
#endif
#ifdef COMPETITIVE_SCHEDULE
  s = new Scheduler<EnqFlags::UPDATEABLE, true>(n);
#else
  s = new Scheduler<EnqFlags::MAYSPEC, true>(n);
#endif
  Update<vertex> u = {GA};
  enqueue_all<EnqFlags(NOHINT | ONLY_ENQUEUES)>(
          swarm::u64it(0), swarm::u64it(n), [&u] (Timestamp, uint64_t v){
          enqueue(init<vertex>, 0ul, s->hint(v), &u, v);
          }, 0ul);

    // Queue each Vertex, prioritized by its current (initial) degree/cardinality
    // FIXME(mcj) as always, we really should invest in a swarm::sort
    std::vector<uint64_t> sortedVertices;
    sortedVertices.reserve(GA.n);
    for (size_t s = 0; s < GA.n; s++) {
        uintE c = GA.V[s].getOutDegree();
        if (c) sortedVertices.push_back((c << 32) | s);
    }
    std::sort(sortedVertices.begin(), sortedVertices.end());


  enqueue_all_progressive<swarm::max_children>(
          sortedVertices.begin(), sortedVertices.end(), [&] (Timestamp ts, uint64_t cv){
                uintE v = cv & ((1ul << 32) - 1);
#ifndef ACCUMULATE
                if (s->extract_ts(v) == u.GA.V[v].getOutDegree())
                    relative_enqueue(s, callUpdate<vertex>, 0ul, v, &u); }, 
#else
                swarm::enqueue(callUpdate<vertex>, cv >> 32, {v, EnqFlags::UPDATEABLE}, &u); },
#endif
          [] (uint64_t cv) { return cv >> 32; }, 
#ifndef ACCUMULATE
          [] (uint64_t v) { return EnqFlags::NOHINT; }); 
#else
          [] (uint64_t v) { return EnqFlags(NOHINT | ONLY_ENQUEUES); });
#endif
  swarm::run();
  s->extract_all_ts(D.s);

  delete s;
  return D;
}

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  cout << "### Application: k-core" << endl;
  cout << "### Graph: " << P.getArgument(0) << endl;
  cout << "### n: " << GA.n << endl;
  cout << "### m: " << GA.m << endl;

  auto cores = KCore(GA);
  uintE mc = 0;
  for (size_t i=0; i < GA.n; i++) { mc = std::max(mc, cores[i]); }
  cout << "### Max core: " << mc << endl;
}
