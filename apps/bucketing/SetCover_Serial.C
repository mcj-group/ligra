// [mcj] Use the AdjacencyGraph format (e.g. in inputs/pbbs/mis/*)

#include "ligra.h"
#include "index_map.h"

#include "swarm/hooks.h"

#include <boost/heap/fibonacci_heap.hpp>
#include <limits>
#include <vector>

// Inexplicably, this is how we make a max heap with the Boost fibonacci_heap
// rather than std::greater. I don't get it...
using PQElement = std::tuple<intE, uintE>;
using PQ = boost::heap::fibonacci_heap<
            PQElement,
            boost::heap::compare<std::less<PQElement> >
            >;

static constexpr uintE THRESHOLD = 2;

template <class vertex>
void SetCover(graph<vertex>& G) {
    static_assert(THRESHOLD >= 1, "We assume the degree threshold is positive");

    timer t; t.start();

    // TODO(mcj):
    // 1) Verify what is the right spot for the zsim_roi_begin and end?

    std::vector<uintE> coveringSet;

    // Conceptually, Approx Set Cover repeatedly removes edges from vertices.
    // We don't actually want to destroy the graph. Originally we used an
    // auxiliary "current degree" vector. But since the algorithm wants to
    // process each vertex in decreasing "current degree" order, let's
    // (1) use a mergeable heap to priority vertices
    // (2) exploit the priority (vertex degree) within the heap to provide the
    // current degree, avoiding extra state.
    PQ pq;
    std::vector<PQ::handle_type> pqHandles;

    coveringSet.reserve(G.n);
    pqHandles.resize(G.n);

    zsim_roi_begin();


    // Queue each vertex, prioritized by its current (initial) degree
    for (size_t v = 0; v < G.n; v++) {
        intE degree = G.V[v].getOutDegree();
        auto handle = pq.push(std::make_tuple(degree, v));
        pqHandles[v] = handle;
    }

    while (!pq.empty()) {
        // pop max-degree vertex
        intE md;
        uintE mdv;
        std::tie(md, mdv) = pq.top();

        if (md <= THRESHOLD) break; // TODO(mcj) make the threshold configurable

        pq.decrease(pqHandles[mdv], std::make_tuple(-1, mdv));
        coveringSet.push_back(mdv);

        // delete set members from universe
        size_t d = G.V[mdv].getOutDegree();
        for (size_t av=0; av < d; av++) {
            uintE adj = G.V[mdv].getOutNeighbor(av);
            auto& handle = pqHandles[adj];
            intE degree = std::get<0>(*handle);
            assert(adj == std::get<1>(*handle));
            if (degree >= 0) {
                intE newDegree = (degree > THRESHOLD) ? degree - 1 : -1;
                pq.decrease(handle, std::make_tuple(newDegree, adj));
            }
        }
    }
    zsim_roi_end();

    t.stop(); t.reportTotal("Running time: ");
    cout << "the covering set length is: " << coveringSet.size() << endl;
//    cout << "the set: ";
//    for (size_t i; i<coveringSetLen; i++) cout << coveringSet[i] << ", ";
    cout << endl;
    return;
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
