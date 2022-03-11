#include "ligra.h"
#include "index_map.h"
#include "bucket.h"
#include "edgeMapReduce.h"

#include "swarm/hooks.h"

template <class vertex>
void SetCover(graph<vertex>& G, size_t num_buckets=128) {
    timer t; t.start();


    size_t n = G.n, m = G.m;
    auto D = array_imap<uintE>(
            G.n,
            [&](size_t i) { return G.V[i].getOutDegree(); });

    // TODO(mcj):
    // 1) What is the difference between Ds below and D above? Aren't they both
    // arrays with the same content?
    // 2) Add a comment along the lines that the reason we add an auxiliary
    // data structure for degrees, rather than just use G.V[i].getOutDegree(),
    // is that the alorithm needs to reduce the degree value, but we don't
    // necessarily want to destroy the graph G itself. (Is that the idea?)
    // 3) Verify what is the right spot for the zsim_roi_begin and end?

    zsim_roi_begin();

    // find degrees
    int Ds[n] = {0}, coveringSet[n], coveringSetLen = 0;
    for (size_t v; v<n; v++) {
        Ds[v] = D(v);
    }

    while (true) {
        // find biggest degree
        int md = -1, mdv=0;
        for (size_t v=0; v<n; v++) {
            if (Ds[v] > md) { md = Ds[v]; mdv = v; }
        }

        if (md <= 2) break; // should add threshold here

        else {
            coveringSet[coveringSetLen++] = mdv;
            Ds[mdv] = -1;
            // delete set members from universe
            for (size_t av=0; av<D(mdv); av++) {
                uintE adj = G.V[mdv].getOutNeighbor(av);
                if (Ds[adj] < 0) continue;
                else if (Ds[adj] > 1) Ds[adj]--;
                else Ds[adj] = -1;
            }
        }
    }
    zsim_roi_end();

    t.stop(); t.reportTotal("Running time: ");
    cout << "the covering set length is: " << coveringSetLen << endl;
//    cout << "the set: ";
//    for (size_t i; i<coveringSetLen; i++) cout << coveringSet[i] << ", ";
    cout << endl;
    return;
}

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
    size_t num_buckets = P.getOptionLongValue("-nb", 128);
    cout << "### Application: set-cover" << endl;
    cout << "### Graph: " << P.getArgument(0) << endl;
    cout << "### Workers: " << getWorkers() << endl;
    cout << "### Buckets: " << num_buckets << endl;
    cout << "### n: " << GA.n << endl;
    cout << "### m: " << GA.m << endl;
    SetCover(GA, num_buckets);
}
