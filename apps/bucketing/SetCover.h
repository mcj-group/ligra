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

    // Validate feasibility of the solution: are all elements at least covered?
    // N.B. validating optimality is hard (this is NP comlete). Maybe we could
    // check the H_n-approximate guarantee? But then we'd need to know the
    // optimal solution...
    for (auto s : cover) {
        size_t sD = G.V[s].getOutDegree();
        for (size_t i = 0; i < sD; i++) {
            isElemCovered[G.V[s].getOutNeighbor(i)] = true;
        }
    }

    auto isTrue = [] (bool b) { return b; };
    bool ec = std::all_of(isElemCovered.begin(), isElemCovered.end(), isTrue);
    bool size = cover.size() <= G.n;
    if (!ec) std::cerr << "ERROR: not all elements were covered" << std::endl;
    if (!size) std::cerr << "ERROR: the subcover was too large" << std::endl;
    if (ec && size) std::cout << "Validation succeeded" << std::endl;
    return ec && size;
}

}
