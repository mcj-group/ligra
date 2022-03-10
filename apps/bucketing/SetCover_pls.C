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
static swarm::bitset isSetCovered;
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


template <class vertex>
static inline void init(swarm::Timestamp, uintE s) {
    uintE c = V<vertex>(s).getOutDegree();
    cardinalities[s] = c;
    if (!c) {
        // Set Cover is not really meant to handle elements not in a set,
        // so handle them semi gracefully by just claiming they are covered.
        // A null set is valid in a set cover, but also uninteresting?
        // I don't think it needs to be added to the subcover
        isSetCovered.set(s, true);
        isElemCovered.set(s, true);
    }
}


template <class vertex>
static inline void coverElement(swarm::Timestamp ts, uintE s, uint64_t elem);


template <class vertex>
static inline void removeSet(swarm::Timestamp ts, uintE s);


template <class vertex>
static inline void addSet(swarm::Timestamp ts, uintE s) {

    // TODO(mcj) could we make the addSet enqueue the correctly timestamped
    // version?
    if (cardinality(ts) > cardinalities[s]) return;

    DEBUG("Add s=%u |s|=%u Deg(s)=%u to the cover\n",
          s, cardinalities[s], V<vertex>(s).getOutDegree());

    // FIXME(mcj) should the cover data structure just be a bit vector that
    // we then collapse at the end?
    cover[s] = s;
    isSetCovered.set(s, true);

    // Delete Set v's member Elements from other Sets
    size_t sD = V<vertex>(s).getOutDegree();

    swarm::enqueue_all_progressive<swarm::max_children>(
        swarm::u64it(0),
        swarm::u64it(sD),
        [s] (swarm::Timestamp ts, uint64_t i) {
             uintE elem = V<vertex>(s).getOutNeighbor(i);
             swarm::enqueue(coverElement<vertex>, ts, EnqFlags::NOHINT, s, elem);
        },
        [ts] (uint64_t) { return ts; },
        [] (uint64_t) { return EnqFlags::NOHINT; }
    );
}


template <class vertex>
static inline void coverElement(swarm::Timestamp ts, uintE s, uintE elem) {
    
    if (isElemCovered.test(elem)) return;

    isElemCovered.set(elem, true);
    size_t elemD = V<vertex>(elem).getOutDegree();

    swarm::enqueue_all_progressive<swarm::max_children>(
        swarm::u64it(0),
        swarm::u64it(elemD),
        [s,elem] (swarm::Timestamp ts, uintE j) {
             uintE s1 = V<vertex>(elem).getOutNeighbor(j);
             if (s1 != s) {
                 swarm::enqueue(removeSet<vertex>, ts, EnqFlags::NOHINT, s1);
             }
         },
        [ts] (uint64_t) { return ts; },
        [] (uint64_t) { return EnqFlags::NOHINT; }
    );
}


template <class vertex>
static inline void removeSet(swarm::Timestamp ts, uintE s) {

    assert(!isSetCovered.test(s));
    uintE card = cardinalities[s];
    assert(card);
    cardinalities[s] -= 1;
    if (card > 1) {
        swarm::Timestamp ts = timestamp(cardinalities[s], s);
        swarm::enqueue(addSet<vertex>, ts, EnqFlags::NOHINT, s);
    } else {
        isSetCovered.set(s, true);
    }
}


template <class vertex>
void SetCover(graph<vertex>& G) {
    vertices = G.V;
    cover.resize(G.n);
    cardinalities.resize(G.n);

    swarm::fill(cover.begin(), cover.end(), INVALID, 0ul);
    isSetCovered.resize<>(G.n, false, 0ul);
    isElemCovered.resize<>(G.n, false, 0ul);
    // Initialize the cardinalities array
    swarm::enqueue_all_progressive<swarm::max_children>(
        swarm::u64it(0),
        swarm::u64it(G.n),
        [] (swarm::Timestamp ts, uint64_t s) {
            swarm::enqueue(init<vertex>, ts,
                           // [mcj] Technically we should use the bit vector
                           // for a hint, as it has a smaller hint space,
                           // but I expect zero-degree vertices to be rare,
                           // so let's not bother
                           swarm::Hint::cacheLine(&cardinalities[s]), s);
        },
        [] (uint64_t) { return 1ul; },
        [] (uint64_t s) { return swarm::Hint::cacheLine(&cardinalities[s]); }
    );

    // Queue each Set, prioritized by its current (initial) degree/cardinality
    // FIXME(mcj) as always, we really should invest in a swarm::sort
    using PQElement = std::tuple<swarm::Timestamp, uintE>;
    std::vector<PQElement> sortedSets;
    sortedSets.reserve(G.n);
    for (size_t s = 0; s < G.n; s++) {
        uintE c = G.V[s].getOutDegree();
        if (c) sortedSets.push_back(std::make_tuple(timestamp(c, s), s));
    }
    std::sort(sortedSets.begin(), sortedSets.end());

    swarm::enqueue_all_progressive<swarm::max_children>(
        sortedSets.begin(),
        sortedSets.end(),
        [] (PQElement& elem) {
            swarm::Timestamp ts;
            uintE s;
            std::tie(ts, s) = elem;
            swarm::enqueue(addSet<vertex>, ts, EnqFlags::NOHINT, s);
        },
        [] (PQElement& elem) { return std::get<0>(elem); },
        [] (PQElement&) -> swarm::Hint { return EnqFlags::NOHINT; }
    );

    swarm::run();

    // FIXME(mcj) add the sort to the ROI. But this requires a parallel Swarm
    // sort, which we haven't yet sorted out.
    std::sort(cover.begin(), cover.end());
    size_t invalids = std::count(cover.begin(), cover.end(), INVALID);
    cover.resize(cover.size() - invalids);

    cout << "the covering set length is: " << cover.size() << endl;

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
