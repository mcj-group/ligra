// [mcj] Use the AdjacencyGraph format (e.g. in inputs/pbbs/mis/*)

#include "ligra.h"
#include "SetCover.h"
#include "swarm/api.h"
#include "swarm/algorithm.h"
#include "swarm/bitset.h"
#include "swarm/scheduler.h"

#include <cassert>
#include <cstdio>
#include <limits>
#include <vector>
#ifdef NONATOMIC_TASKS
#include <atomic>
#endif

#define DEBUG(args...) //swarm::info(args)

static constexpr uintE INVALID = std::numeric_limits<uintE>::max();
#ifdef COMPETITIVE_SCHEDULE
//static swarm::Scheduler<EnqFlags::UPDATEABLE, true> *sets;
static swarm::Scheduler<EnqFlags::COMPETITIVE, false> *elements;
#elif defined(HIVE_BASIC)
static swarm::Scheduler<EnqFlags::UPDATEABLE, true> *elements;
#endif
//#else
#ifdef NONATOMIC_TASKS
static std::vector<std::atomic<int64_t>>* cardinalities;
static std::vector<std::atomic_flag>* isElemCovered;
#else
static std::vector<uintE> cardinalities;
static swarm::bitset isElemCovered;
#endif
//#endif

static std::vector<uintE> cover;
// [mcj] Major hack to get around the annoying templating of the vertex type
static void* vertices;

template <class vertex>
const vertex& V(uint64_t v) { return reinterpret_cast<vertex*>(vertices)[v]; }


constexpr swarm::Timestamp timestamp(uint64_t cardinality, uint64_t setID) {
    // Process sets in decreasing order of cardinality
    uint64_t invert = UINT32_MAX - 8ul - cardinality;
    return (invert << 32ul) | setID;
}


constexpr uintE cardinality(swarm::Timestamp ts) {
    return UINT32_MAX - 8ul - (ts >> 32);
}


constexpr uintE set(swarm::Timestamp ts) {
    return ts & ((1ul << 32) - 1);
}


static inline uint64_t hint(uintE set) {
#if defined(COMPETITIVE_SCHEDULE) || defined(HIVE_BASIC)
    return set;
#elif defined(NONATOMIC_TASKS)
    return swarm::Hint::cacheLine(&(*cardinalities)[set]);
#else
    return swarm::Hint::cacheLine(&cardinalities[set]);
#endif
}


template <class vertex>
static inline void init(swarm::Timestamp, uintE s) {
//#ifdef COMPETITIVE_SCHEDULE
//  sets->set_base_ts(s, timestamp(V<vertex>(s).getOutDegree(), s));
//#else
#ifdef NONATOMIC_TASKS
    std::atomic_init(&(*cardinalities)[s], V<vertex>(s).getOutDegree());
#else
    cardinalities[s] = V<vertex>(s).getOutDegree();
#endif
//#endif
}

template <class vertex>
static inline void coverElement(swarm::Timestamp ts, uintE s, uint64_t elem);

template <class vertex>
static inline void addSet(swarm::Timestamp ts, uintE s) {
//#ifndef COMPETITIVE_SCHEDULE
#ifdef NONATOMIC_TASKS
    int64_t card = (*cardinalities)[s].load(std::memory_order_relaxed);
    if (cardinality(ts) > card) {
        // This task instance is too early given the now-lower |s|
        if (card > 0) {
            // So re-enqueue with the correct timestamp for |s|.
            ts = timestamp(card, s);
#else
    if (cardinality(ts) > cardinalities[s]) {
        // This task instance is too early given the now-lower |s|
        if (cardinalities[s]) {
            // So re-enqueue with the correct timestamp for |s|.
            ts = timestamp(cardinalities[s], s);
#endif
            EnqFlags flags = EnqFlags(SAMEHINT | SAMETASK | MAYSPEC);
            swarm::enqueue(addSet<vertex>, ts, flags, s);
        }
        return;
    }
//#endif

    DEBUG("Add s=%u |s|=%u Deg(s)=%u to the cover\n",
          s, cardinality(ts), V<vertex>(s).getOutDegree());

    // FIXME(mcj) should the cover data structure just be a bit vector that
    // we then collapse at the end?
    cover[s] = s;
    // s and all its elements are now covered, so reset s's effective
    // cardinality to filter away all other tasks for s.
//#ifdef COMPETITIVE_SCHEDULE
//    sets->set_base_ts(s, swarm::NEVER);
//#else
#ifdef NONATOMIC_TASKS
    (*cardinalities)[s].store(0, std::memory_order_relaxed);
#else
    cardinalities[s] = 0;
#endif
//#endif

    // Delete Set v's member Elements from other Sets
    const vertex& vs = V<vertex>(s);
    size_t sD = vs.getOutDegree();
#if 0
    bool redundant = true;
    for (int i = 0; i < sD; i++) {
        if (!isElemCovered.test(vs.getOutNeighbor(i))) {
            redundant = false;
            break;
        }
    }
    if (redundant) {
        swarm::info("added redundant set to cover with covered neighbours:");
        for (int i = 0; i < sD; i++) {
            swarm::info("%u", vs.getOutNeighbor(i));
        }
        assert(false);
    }
#endif
    swarm::enqueue_all<EnqFlags(NOHINT | MAYSPEC)>(
        swarm::u64it(0),
        swarm::u64it(sD),
        [s,&vs] (swarm::Timestamp ts, uint64_t i) {
            uintE elem = vs.getOutNeighbor(i);
#if defined(COMPETITIVE_SCHEDULE)
            swarm::absolute_enqueue(elements, coverElement<vertex>,
                    ts, elem, s, elem);
#elif defined(HIVE_BASIC)
            //DEBUG("prev ts: %lu, curr ts: %lu",
            //        element->extract_ts(elem), ts);
            if (elements->extract_ts(elem) > ts)
                swarm::absolute_enqueue(elements, coverElement<vertex>,
                    ts, elem, s, elem);
#else
            swarm::enqueue(coverElement<vertex>, ts,
#ifdef NONATOMIC_TASKS
                           {swarm::Hint::cacheLine(&((*isElemCovered)[elem])), EnqFlags::MAYSPEC},
#else
                           {isElemCovered.hint(elem), EnqFlags(PRODUCER | MAYSPEC)},
#endif
                           s, elem);
#endif
        },
        ts);
}


/*#ifdef COMPETITIVE_SCHEDULE
template<class vertex>
static inline void decrementCardinality(swarm::Timestamp, uintE s) {
    if (pls_unlikely(cardinality(sets->extract_ts(s)) == 1)) {
        swarm::absolute_enqueue(sets, addSet<vertex>, swarm::NEVER, int(s), s);
    }
    else {
        swarm::relative_enqueue(sets, addSet<vertex>, 1ul << 32, hint(s), s);
    }
}
#else*/
template<class vertex>
#ifdef NONATOMIC_TASKS
static inline void decrementCardinality(swarm::Timestamp, std::atomic<int64_t>* cptr) {
    DEBUG("decrement cardinality of set %lu to %lu", std::distance(&cardinalities[0], cptr), *cptr - 1);
    cptr->fetch_sub(1, std::memory_order_relaxed);
}
#else
static inline void decrementCardinality(swarm::Timestamp, uintE* cptr) {
    DEBUG("decrement cardinality of set %lu to %lu", std::distance(&cardinalities[0], cptr), *cptr - 1);
    if(*cptr > 0) (*cptr)--;
}
#endif
#if 0
    if(*cptr == 0) {
        const vertex& s = V<vertex>(std::distance(&cardinalities[0], cptr));
        for (int i = 0; i < s.getOutDegree(); i++) {
            if (!isElemCovered.test(s.getOutNeighbor(i)))
                swarm::info("uncovered neighbor of set %lu with 0 cardinality: %lu",
                        std::distance(&cardinalities[0], cptr), s.getOutNeighbor(i));
            assert(isElemCovered.test(s.getOutNeighbor(i)));
        }
    }
#endif
//}
//#endif


template <class vertex>
static inline void coverElement(swarm::Timestamp ts, uintE s, uintE elem) {
    DEBUG("%lu: Cover element %lu by set %lu, degree %u, currently %scovered",
            ts, elem, s, V<vertex>(elem).getInDegree(), isElemCovered.test(elem) ? "" : "un");
#if !defined(COMPETITIVE_SCHEDULE) && !defined(HIVE_BASIC)
#ifdef NONATOMIC_TASKS
    if ((*isElemCovered)[elem].test_and_set()) return;
#else
    if (isElemCovered.test(elem)) return;
    isElemCovered.set(elem, true);
#endif
#endif

    const vertex& ve = V<vertex>(elem);
    size_t elemD = ve.getInDegree();
    swarm::enqueue_all<EnqFlags(NOHINT | MAYSPEC)>(
        swarm::u64it(0),
        swarm::u64it(elemD),
        [s,&ve] (swarm::Timestamp ts, uintE j) {
            uintE s1 = ve.getInNeighbor(j);
            if (s1 != s) {
/*#ifdef COMPETITIVE_SCHEDULE
                swarm::enqueue(decrementCardinality<vertex>, ts, hint(s1), s1);
#else*/
#ifdef NONATOMIC_TASKS
                auto cptr = &((*cardinalities)[s1]);
#else
                auto cptr = &cardinalities[s1];
#endif
                swarm::enqueue(decrementCardinality<vertex>, ts,
                        {hint(s1), EnqFlags::MAYSPEC}, cptr);
//#endif
            }
        },
        ts);
}


template <class vertex>
static inline void addSetCG(swarm::Timestamp ts, uintE s) {
#ifndef NONATOMIC_TASKS
#if !defined(COMPETITIVE_SCHEDULE) && !defined(HIVE_BASIC)
    // FIXME the CG version doesn't need to break ties for equal-cardinality,
    // so this is unnecessarily deterministic
    if (cardinality(ts) > cardinalities[s]) {
        // This task instance is too early given the now-lower |s|
        if (cardinalities[s]) {
            // So re-enqueue with the correct timestamp for |s|.
            ts = timestamp(cardinalities[s], s);
            EnqFlags flags = EnqFlags(SAMEHINT | SAMETASK);
            swarm::enqueue(addSetCG<vertex>, ts, flags, s);
        }
        return;
    }
#endif

    DEBUG("Add s=%u |s|=%u Deg(s)=%u to the cover\n",
          s, cardinality(ts), V<vertex>(s).getOutDegree());

    // FIXME(mcj) should the cover data structure just be a bit vector that
    // we then collapse at the end?
    cover[s] = s;
    // s and all its elements are now covered, so reset s's effective
    // cardinality to filter away all other tasks for s.
//#ifdef COMPETITIVE_SCHEDULE
//    sets->set_base_ts(s, swarm::NEVER);
//#else
    cardinalities[s] = 0;
//#endif

    // Delete Set v's member Elements from other Sets
    const vertex& vs = V<vertex>(s);
    size_t sD = vs.getOutDegree();
    for (size_t i = 0; i < sD; i++) {
        uintE elem = vs.getOutNeighbor(i);
        if (isElemCovered.test(elem)) continue;
        isElemCovered.set(elem, true);
        const vertex& ve = V<vertex>(elem);
        size_t elemD = ve.getInDegree();
        for (size_t j = 0; j < elemD; j++) {
            uintE s1 = ve.getInNeighbor(j);
            if (s1 != s) {
//#ifdef COMPETITIVE_SCHEDULE
                //TODO: Make this an enqueuer tree to cg_comp will actually work
//                swarm::relative_enqueue(sets, addSetCG<vertex>, 1ul << 32, hint(s1), s1);
//#else
                cardinalities[s1]--;
//#endif
            }
        }
    }
#endif
}


template <class vertex>
void SetCover(graph<vertex>& G) {
    vertices = G.V;
    cover.resize(G.n);
#ifdef COMPETITIVE_SCHEDULE
  //  sets = new swarm::Scheduler<EnqFlags::UPDATEABLE, true>(G.n);
    elements = new swarm::Scheduler<EnqFlags::COMPETITIVE, false>(G.n);
#elif defined(HIVE_BASIC)
    elements = new swarm::Scheduler<EnqFlags::UPDATEABLE, true>(G.n);
#endif
//#else
#ifndef NONATOMIC_TASKS
    cardinalities.resize(G.n);
#endif
#ifdef RELAXED
    //These need to happen before everything else, which
    //I can't guarantee with the relaxed scheduler
    //so they need to be moved before the ROI
    std::fill(cover.begin(), cover.end(), INVALID);
#ifdef NONATOMIC_TASKS
    std::vector<std::atomic<int64_t>> cards(G.n);
    cardinalities = &cards;
    std::vector<std::atomic_flag> elems(G.n);
    for (int i = 0; i < G.n; i++) elems[i].clear();
    isElemCovered = &elems;
#else
    isElemCovered.resize(G.n, false);
#endif
    for (int i = 0; i < G.n; i++) {
        init<vertex>(0, i);
    }
#else
//#endif
    swarm::fill(cover.begin(), cover.end(), INVALID, 0ul);
    isElemCovered.resize<>(G.n, false, 0ul);
    // Initialize the cardinalities array
    swarm::enqueue_all_progressive<swarm::max_children>(
        swarm::u64it(0),
        swarm::u64it(G.n),
        [] (swarm::Timestamp ts, uint64_t s) {
            swarm::enqueue(init<vertex>, ts, hint(s), s);
        },
        [] (uint64_t) { return 1ul; },
        [] (uint64_t s) { return hint(s); }
    );
#endif

    // Queue each Set, prioritized by its current (initial) degree/cardinality
    // FIXME(mcj) as always, we really should invest in a swarm::sort
    std::vector<swarm::Timestamp> sortedSets;
    sortedSets.reserve(G.n);
    for (size_t s = 0; s < G.n; s++) {
        uintE c = G.V[s].getOutDegree();
        if (c) sortedSets.push_back(timestamp(c, s));
    }
    std::sort(sortedSets.begin(), sortedSets.end());

    swarm::enqueue_all_progressive<swarm::max_children>(
        sortedSets.begin(),
        sortedSets.end(),
        [] (swarm::Timestamp ts) {
            uintE s = set(ts);
//#ifdef COMPETITIVE_SCHEDULE
            //swarm::info("set = %u, ts = %lu, setTS = %lu", s, ts, sets->extract_ts(s));
//            swarm::relative_enqueue(sets,
//#else
            swarm::enqueue(
//#endif
#ifdef COARSE_GRAIN
                           addSetCG<vertex>,
#else
                           addSet<vertex>,
#endif
//#ifdef COMPETITIVE_SCHEDULE
//                          0,
//#else
                           ts,
//#endif
                           {hint(s), EnqFlags::MAYSPEC}, s);
        },
        [] (swarm::Timestamp ts) { return ts; },
        [] (swarm::Timestamp) -> swarm::Hint { return EnqFlags(NOHINT | MAYSPEC); }
    );
    void* ifn = reinterpret_cast<void*>(swarm::bareRunner<decltype(init<vertex>), init<vertex>, uintE>);
    swarm::programTSP(ifn, 10, 0, reinterpret_cast<uintptr_t>(vertices), 32, 8);
    swarm::programTSP(ifn, 11, 0, reinterpret_cast<uintptr_t>(&cardinalities[0]), 8, 8);
    void* afn = reinterpret_cast<void*>(swarm::bareRunner<decltype(addSet<vertex>), addSet<vertex>, uintE>);
    swarm::programTSP(afn, 10, 0, reinterpret_cast<uintptr_t>(&cardinalities[0]), 8, 8);
    swarm::programTSP(afn, 11, 0, reinterpret_cast<uintptr_t>(&cover[0]), 8, 8);
    swarm::programTSP(afn, 12, 0, reinterpret_cast<uintptr_t>(vertices), 32, 8);
    void* cfn = reinterpret_cast<void*>(swarm::bareRunner<decltype(coverElement<vertex>), coverElement<vertex>, uintE, uintE>);
    swarm::programTSP(cfn, 10, 1, reinterpret_cast<uintptr_t>(vertices), 32, 8);
#ifdef NONATOMIC_TASKS
    void* dfn = reinterpret_cast<void*>(swarm::bareRunner<decltype(decrementCardinality<vertex>), decrementCardinality<vertex>, std::atomic<int64_t>*>);
#else
    void* dfn = reinterpret_cast<void*>(swarm::bareRunner<decltype(decrementCardinality<vertex>), decrementCardinality<vertex>, uintE*>);
#endif
    swarm::programTSP(dfn, 10, 0, 0, 1, 8);

    swarm::run();

    // FIXME(mcj) add the sort to the ROI. But this requires a parallel Swarm
    // sort, which we haven't yet sorted out.
    std::sort(cover.begin(), cover.end());
    size_t invalids = std::count(cover.begin(), cover.end(), INVALID);
    cover.resize(cover.size() - invalids);

    if (!setcover::success<vertex>(G, cover)) std::abort();
}


template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
    cout << "### Application: set-cover" << endl;
    cout << "### Graph: " << P.getArgument(0) << endl;
    cout << "### Workers: " << getWorkers() << endl;
    cout << "### n: " << GA.n << endl;
    cout << "### m: " << GA.m << endl;
    SetCover(GA);
}
