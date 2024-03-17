#include "ligra.h"
#include "index_map.h"
#include "bucket.h"
#include "edgeMapReduce.h"
#include "swarm/api.h"
#include "swarm/scheduler.h"
using namespace swarm;

#ifdef COMPETITIVE_SCHEDULE
Scheduler<EnqFlags::COMPETITIVE, true> *s;
#else
Scheduler<EnqFlags::MAYSPEC, true> *s;
#endif

template<class vertex> struct Update;

//function because functors don't seem to work properly with enqueue macros
template<class vertex>
static void callUpdate(Timestamp ts, Update<vertex> *u, uintE v)
{
    (*u)(ts, v);
}

//function to call the relative update for fine-grained implementations
template<class vertex>
static void decrementDegree(Timestamp ts, Update<vertex> *u, uintE v)
{
#ifndef COMPETITIVE_SCHEDULE
    if (s->extract_ts(v) > ts) //early exit simulates COMPETITIVE semantics w/out the flag
#endif
        relative_enqueue(s, callUpdate<vertex>, -1, v, u, v);
}

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
        enqueue_all_progressive<swarm::max_children>(
                swarm::u64it(0), swarm::u64it(GA.V[v].getOutDegree()),
                [this, v] (Timestamp ts, uint64_t i) {
                    uintE ngh = GA.V[v].getOutNeighbor(i);
#ifdef COARSE_GRAIN
#ifndef COMPETITIVE_SCHEDULE
                    decrementDegree<vertex>(ts, this, ngh); },
#else
                    relative_enqueue(s, callUpdate<vertex>, -1, ngh, this, ngh); },
#endif
#else
                    enqueue(decrementDegree<vertex>, ts, {ngh, EnqFlags::MAYSPEC}, this, ngh); },
#endif
                [ts] (uint64_t) { return ts; },
                [] (uint64_t) { return EnqFlags(NOHINT | MAYSPEC); }); //no point in a hint unless graph is well-ordered/low degree
    }
};

template <class vertex>
array_imap<uintE> KCore(graph<vertex>& GA, size_t num_buckets=128) {
  const size_t n = GA.n; const size_t m = GA.m;
  auto D = array_imap<uintE>(n, [&] (size_t i) { return GA.V[i].getOutDegree(); });
#ifdef COMPETITIVE_SCHEDULE
  s = new Scheduler<EnqFlags::COMPETITIVE, true>(n);
#else
  s = new Scheduler<EnqFlags::MAYSPEC, true>(n);
#endif
  Update<vertex> u = {GA};
  enqueue_all_progressive<swarm::max_children>(
          swarm::u64it(0), swarm::u64it(n), [&] (Timestamp ts, uint64_t v){
                absolute_enqueue(s, callUpdate<vertex>, GA.V[v].getOutDegree(), v, &u, v); },
          [] (uint64_t) { return 0ul; },
          [] (uint64_t v) { return swarm::Hint({v, EnqFlags::MAYSPEC}); }); //accessing vertices in order so might as well use as hint
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
