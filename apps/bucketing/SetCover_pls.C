// [mcj] Use the AdjacencyGraph format (e.g. in inputs/pbbs/mis/*)

#include "ligra.h"
#include "SetCover.h"
#include "swarm/api.h"
#include "swarm/algorithm.h"
#include "swarm/bitset.h"

#include <cassert>
#include <cstdio>
#include <limits>
#include <vector>

#define DEBUG(args...) //swarm::info(args)

static constexpr uintE INVALID = std::numeric_limits<uintE>::max();


static std::vector<uintE> cardinalities;
static std::vector<uintE> cover;
static swarm::bitset isElemCovered;
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
    return swarm::Hint::cacheLine(&cardinalities[set]);
}


template <class vertex>
static inline void init(swarm::Timestamp, uintE s) {
    cardinalities[s] = V<vertex>(s).getOutDegree();
}


template <class vertex>
static inline void coverElement(swarm::Timestamp ts, uintE s, uint64_t elem);


template <class vertex>
static inline void addSet(swarm::Timestamp ts, uintE s) {
    if (cardinality(ts) > cardinalities[s]) {
        // This task instance is too early given the now-lower |s|
        if (cardinalities[s]) {
            // So re-enqueue with the correct timestamp for |s|.
            ts = timestamp(cardinalities[s], s);
            EnqFlags flags = EnqFlags(SAMEHINT | SAMETASK);
            swarm::enqueue(addSet<vertex>, ts, flags, s);
        }
        return;
    }

    DEBUG("Add s=%u |s|=%u Deg(s)=%u to the cover\n",
          s, cardinalities[s], V<vertex>(s).getOutDegree());

    // FIXME(mcj) should the cover data structure just be a bit vector that
    // we then collapse at the end?
    cover[s] = s;
    // s and all its elements are now covered, so reset s's effective
    // cardinality to filter away all other tasks for s.
    cardinalities[s] = 0;

    // Delete Set v's member Elements from other Sets
    const vertex& vs = V<vertex>(s);
    size_t sD = vs.getOutDegree();
    swarm::enqueue_all<EnqFlags::NOHINT>(
        swarm::u64it(0),
        swarm::u64it(sD),
        [s,&vs] (swarm::Timestamp ts, uint64_t i) {
            uintE elem = vs.getOutNeighbor(i);
            swarm::enqueue(coverElement<vertex>, ts,
                           {isElemCovered.hint(elem), EnqFlags::PRODUCER},
                           s, elem);
        },
        ts);
}


template <class vertex>
static inline void decrementCardinality(swarm::Timestamp, uintE* cptr) {
    (*cptr)--;
}


template <class vertex>
static inline void coverElement(swarm::Timestamp ts, uintE s, uintE elem) {
    if (isElemCovered.test(elem)) return;
    isElemCovered.set(elem, true);

    const vertex& ve = V<vertex>(elem);
    size_t elemD = ve.getInDegree();
    swarm::enqueue_all<EnqFlags::NOHINT>(
        swarm::u64it(0),
        swarm::u64it(elemD),
        [s,&ve] (swarm::Timestamp ts, uintE j) {
            uintE s1 = ve.getInNeighbor(j);
            if (s1 != s) {
                uintE* cptr = &cardinalities[s1];
                swarm::enqueue(decrementCardinality<vertex>, ts, hint(s1), cptr);
            }
        },
        ts);
}


template <class vertex>
static inline void addSetCG(swarm::Timestamp ts, uintE s) {
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

    DEBUG("Add s=%u |s|=%u Deg(s)=%u to the cover\n",
          s, cardinalities[s], V<vertex>(s).getOutDegree());

    // FIXME(mcj) should the cover data structure just be a bit vector that
    // we then collapse at the end?
    cover[s] = s;
    // s and all its elements are now covered, so reset s's effective
    // cardinality to filter away all other tasks for s.
    cardinalities[s] = 0;

    // Delete Set v's member Elements from other Sets
    const vertex& vs = V<vertex>(s);
    size_t sD = vs.getOutDegree();
    for (size_t i = 0; i < sD; i++) {
        uintE elem = vs.getOutNeighbor(i);
        if (isElemCovered.test(elem)) continue;
        isElemCovered.set(elem, true);

        const vertex& ve = V<vertex>(elem);
        size_t elemD = ve.getOutDegree();
        for (size_t j = 0; j < elemD; j++) {
            uintE s1 = ve.getOutNeighbor(j);
            if (s1 != s) {
                cardinalities[s1]--;
            }
        }
    }
}


template <class vertex>
void SetCover(graph<vertex>& G) {
    vertices = G.V;
    cover.resize(G.n);
    cardinalities.resize(G.n);

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
            swarm::enqueue(
#ifdef COARSE_GRAIN
                           addSetCG<vertex>,
#else
                           addSet<vertex>,
#endif
                           ts, hint(s), s);
        },
        [] (swarm::Timestamp ts) { return ts; },
        [] (swarm::Timestamp) -> swarm::Hint { return EnqFlags::NOHINT; }
    );

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
