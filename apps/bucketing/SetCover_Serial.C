// [mcj] Use the AdjacencyGraph format (e.g. in inputs/pbbs/mis/*)

#include "ligra.h"
#include "index_map.h"

#include "swarm/hooks.h"

#include <algorithm>
#include <limits>
#include <vector>

template <class vertex>
void SetCover(graph<vertex>& G) {
    timer t; t.start();

    // TODO(mcj):
    // 1) Verify what is the right spot for the zsim_roi_begin and end?

    // find degrees
    std::vector<uintE> coveringSet;
    coveringSet.reserve(G.n);

    zsim_roi_begin();

    // Approx Set Cover repeatedly removes edges from vertices. We don't
    // actually want to destroy the graph, so we use the following
    // temporary/auxiliary structure to hold the "current" degree of each vertex
    auto D = array_imap<intE>(
            G.n,
            [&](size_t i) -> intE { return G.V[i].getOutDegree(); });

    intE md = std::numeric_limits<intE>::max();
    while (md > 2) { // should add threshold here
        // find biggest degree
        intE* mdp = std::max_element(D.s, D.e);
        md = *mdp;
        uintE mdv = std::distance(D.s, mdp);

        coveringSet.push_back(mdv);
        *mdp = -1;
        D[mdv] = -1;
        // delete set members from universe
        size_t d = G.V[mdv].getOutDegree();
        for (size_t av=0; av < d; av++) {
            uintE adj = G.V[mdv].getOutNeighbor(av);
            if (D[adj] < 0) continue;
            else if (D[adj] > 1) D[adj]--;
            else D[adj] = -1;
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
