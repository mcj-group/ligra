#pragma once

#include "ligra.h"

#include <algorithm>
#include <iostream>
#include <vector>

namespace setcover {

template <class vertex, class collection>
inline bool success(const graph<vertex>& G, const collection& cover) {
    std::vector<bool> isElemCovered(G.n, false);
    // Handle disconnected graphs, just in case
    for (size_t elem = 0; elem < G.n; elem++) {
        if (!G.V[elem].getOutDegree()) isElemCovered[elem] = true;
    }

    // Validate feasibility of the solution:
    // 1) All elements should be covered by the set cover
    // 2) No set in the cover should be redundant (i.e. all its elements are
    //    already covered)
    // 3) N.B. validating optimality is hard (this is NP comlete). Maybe we
    //    could check the H_n-approximate guarantee? But then we'd need to know
    //    the optimal solution...
    std::vector<uint64_t> redundantSets;
    for (auto s : cover) {
        size_t sD = G.V[s].getOutDegree();
        bool redundant = true;
        for (size_t i = 0; i < sD; i++) {
            size_t elem = G.V[s].getOutNeighbor(i);
            if (!isElemCovered[elem]) {
                redundant = false;
                isElemCovered[elem] = true;
            }
        }
        if (redundant) redundantSets.push_back(s);
    }

    auto isTrue = [] (bool b) { return b; };
    bool ec = std::all_of(isElemCovered.begin(), isElemCovered.end(), isTrue);
    bool size = cover.size() <= G.n;
    bool red = !redundantSets.empty();
    if (!ec) std::cerr << "ERROR: not all elements were covered" << std::endl;
    if (!size) std::cerr << "ERROR: the subcover was too large" << std::endl;
    if (red) {
        std::cerr << "ERROR: found redundant sets in the cover: " << std::endl;
        for (auto s : redundantSets) std::cerr << s << std::endl;
    }
    bool success = (ec && size && !red);
    if (success) std::cout << "Validation succeeded" << std::endl;
    return success;
}

}
