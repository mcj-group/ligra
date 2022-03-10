// [mcj] Use the AdjacencyGraph format (e.g. in inputs/pbbs/mis/*)

#include "ligra.h"
#include "SetCover.h"

#include "swarm/hooks.h"

#include <boost/heap/fibonacci_heap.hpp>
#include <cassert>
#include <cstdio>
#include <limits>
#include <vector>

#define DEBUG(args...) //printf(args)

// Inexplicably, this is how we make a max heap with the Boost fibonacci_heap
// rather than std::greater. I don't get it...
using PQElement = std::tuple<uintE, uintE>;
using PQ = boost::heap::fibonacci_heap<
            PQElement,
            boost::heap::compare<std::less<PQElement> >
            >;


template <class vertex>
void SetCover(graph<vertex>& G) {
    timer t; t.start();

    // [mcj] As mentioned in the in Julienne paper, we will generate "bipartite
    // graphs to use as set cover instances by having vertices represent both
    // the sets and the elements." I *think* they meant transforming the graph
    // to a Bipartite Double Cover
    // https://en.wikipedia.org/wiki/Bipartite_double_cover So we have one
    // vector to represent the Element vertices, and another vector to represent
    // the Set vertices. An element is in a given set if there was an edge
    // connecting their corresponding two vertices in the original graph.
    //
    // Johnson's 1974 Approx Set Cover selects the maximum cardinality uncovered
    // set adds it to the cover, and removes its now-covered elements from other
    // sets. Processing Set vertices in decreasing degree order, and removing
    // edges from the Set vertices. But we don't actually want to destroy the
    // original graph.
    //
    // Originally we used an auxiliary "current degree" vector for this reason,
    // but since the algorithm wants to process each vertex in decreasing
    // "current degree" order, let's
    // (1) use a mergeable heap to prioritize Set vertices
    // (2) exploit the priority (Set vertex degree) within the heap to provide the
    // current degree, avoiding extra state.
    //
    // FIXME(mcj) Dhulipala et al. claim Johnson's algorithm has time complexity
    // O(m) for unweighted sets, but by using a Fibonacci heap, the algorithm
    // below looks like O(n lg n + m).
    // Would a Bucket Queue improve the work to linear?

    std::vector<uintE> cover;
    cover.reserve(G.n);
    PQ pq;
    std::vector<PQ::handle_type> pqHandles(G.n);

    zsim_roi_begin();

    std::vector<bool> isSetCovered(G.n, false);
    std::vector<bool> isElemCovered(G.n, false);

    // Queue each Set, prioritized by its current (initial) degree/cardinality
    for (size_t s = 0; s < G.n; s++) {
        intE cardinality = G.V[s].getOutDegree();
        if (cardinality) {
            pqHandles[s] = pq.push({cardinality, s});
        } else {
            // TODO(mcj) Reason more about disconnected graphs.
            // Set Cover is not really meant to handle elements not in a set,
            // so handle them semi gracefully by just claiming they are covered.
            // A null set is valid in a set cover, but also uninteresting?
            // I don't think it needs to be added to the subcover
            isSetCovered[s] = true;
            isElemCovered[s] = true;
        }
    }

    while (!pq.empty()) {
        // Pop the max-cardinality (degree) set (vertex)
        uintE cardinality;
        uintE s;
        std::tie(cardinality, s) = pq.top();
        pq.pop();

        DEBUG("Add s=%u |s|=%u Deg(s)=%u to the cover\n",
              s, cardinality, G.V[s].getOutDegree());
        assert(cardinality);

        cover.push_back(s);
        isSetCovered[s] = true;

        // Delete Set v's member Elements from other Sets
        const vertex& vs = G.V[s];
        size_t sD = vs.getOutDegree();
        for (size_t i = 0; i < sD; i++) {
            uintE elem = vs.getOutNeighbor(i);
            if (isElemCovered[elem]) continue;
            isElemCovered[elem] = true;

            const vertex& ve = G.V[elem];
            size_t elemD = ve.getOutDegree();
            for (size_t j = 0; j < elemD; j++) {
                uintE s1 = ve.getOutNeighbor(j);
                if (s1 == s) continue;

                assert(!isSetCovered[s1]); // Otherwise elem should be covered
                auto& handle = pqHandles[s1];
                uintE card = std::get<0>(*handle);
                assert(s1 == std::get<1>(*handle));
                if (card > 1) {
                    pq.decrease(handle, {card - 1, s1});
                } else {
                    // All of s1's elements have been covered by other sets, so
                    // mark it covered, but do not add it to the subcover
                    pq.erase(handle);
                    isSetCovered[s1] = true;
                }
            }
        }
    }
    zsim_roi_end();

    t.stop(); t.reportTotal("Running time: ");
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
