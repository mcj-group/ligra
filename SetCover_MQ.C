// [mcj] Use the AdjacencyGraph format (e.g. in inputs/pbbs/mis/*)

#include "ligra.h"
#include "SetCover.h"
#include "MultiQueue.h"
#include "MultiBucketQueue.h"

#include <cassert>
#include <cstdio>
#include <thread>
#include <limits>
#include <vector>
#include <functional>
#include <set>

using PQElement = std::tuple<uintE, uintE>;

struct stats {
  uint64_t emptyWork = 0;
};

template <class vertex, typename MQ>
void MQThreadTask(graph<vertex>& G, MQ& wl,
                            atomic<bool>* isElemCovered,
                            atomic<uint32_t>* cardinality,
                            atomic<bool>* cover,
                            stats* threadStat)
{
    uint64_t emptyWork = 0;
    wl.initTID();

    while (true) {
        // Pop the max-cardinality (degree) set (vertex)
        uintE pushedCard, s;
        auto item = wl.pop();
        if (item) std::tie(pushedCard, s) = item.get();
        else break;

        if (pushedCard == 0) {
            emptyWork++;
            continue;
        }

        // postponed and all its elements are already covered
        // should not be added to subcover
        uint32_t curCard = cardinality[s].load(std::memory_order_seq_cst);
        if (curCard == 0) {
            emptyWork++;
            continue;
        }

        // postponed
        if (pushedCard > curCard) {
            wl.push(curCard, s);
            emptyWork++;
            continue;
        }

        // check if this is already in subcover
        bool flag = false;
        if (!cover[s].compare_exchange_weak(flag, true, memory_order_seq_cst))
            continue;
            
        // Delete Set v's member Elements from other Sets
        const vertex& vs = G.V[s];
        size_t sD = vs.getOutDegree();
        for (size_t i = 0; i < sD; i++) {
            uintE elem = vs.getOutNeighbor(i);

            // check if this node's member elements
            // have already been processed
            bool processed = false;
            if (!isElemCovered[elem].compare_exchange_weak(
                    processed, true, memory_order_seq_cst))
                continue;

            const vertex& ve = G.V[elem];
            size_t elemD = ve.getInDegree();
            for (size_t j = 0; j < elemD; j++) {
                uintE s1 = ve.getInNeighbor(j);
                if (s1 == s) continue;

                // decrease the outCardinlaity (priority)
                uint32_t card = cardinality[s1].load(std::memory_order_seq_cst);
                bool decreased = false;
                while (!decreased && card > 0) {
                    decreased = cardinality[s1].compare_exchange_weak(
                        card, card - 1,
                        memory_order_seq_cst, memory_order_seq_cst);
                };
            }
        }
    }

    threadStat->emptyWork = emptyWork;
}

template <class vertex, typename MQ_Type>
void spawnTasks(graph<vertex>& G, MQ_Type &wl, int threadNum, 
                atomic<bool>* isElemCovered, atomic<uint32_t>* cardinality,
                atomic<bool>* cover, bool noverify=false)
{
    int cnt1 = 0, cnt2 = 0;
    // Queue each Set, prioritized by its current (initial) degree/cardinality
    for (size_t s = 0; s < G.n; s++) {
        uintE outCardinality = G.V[s].getOutDegree();
        uintE inCardinality = G.V[s].getInDegree();

        if (outCardinality > 0) {
            wl.push(outCardinality, s);
            cardinality[s] = outCardinality;
        }

        if (outCardinality == 0) cnt1++;
        if (inCardinality == 0) cnt2++;
    }

    cout << cnt1 << " vertices with no outgoing edge\n";
    cout << cnt2 << " vertices with no incoming edge\n";

    stats threadStats[threadNum];
    for (int i = 0; i < threadNum;i++) {
        threadStats[i].emptyWork = 0;
    }

    auto begin = std::chrono::high_resolution_clock::now();

    vector<thread*> workers;
    cpu_set_t cpuset;
    for (int i = 1; i < threadNum; i++) {
        CPU_ZERO(&cpuset);
        uint64_t coreID = i;
        CPU_SET(coreID, &cpuset);
        std::thread *newThread = new std::thread(
            MQThreadTask<vertex, MQ_Type>, ref(G),
            ref(wl), ref(isElemCovered),
            ref(cardinality), ref(cover),
            &threadStats[i]
        );
        int rc = pthread_setaffinity_np(newThread->native_handle(),
                                        sizeof(cpu_set_t), &cpuset);
        if (rc != 0) {
            std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
        }
        workers.push_back(newThread);
    }
    CPU_ZERO(&cpuset);
    CPU_SET(0, &cpuset);
    sched_setaffinity(0, sizeof(cpuset), &cpuset);
    MQThreadTask<vertex, MQ_Type>(G, wl, isElemCovered, cardinality, cover, &threadStats[0]);
    for (thread*& worker : workers) {
        worker->join();
        delete worker;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
    wl.stat();
    std::cout << "runtime_ms " << ms << "\n";

    uint64_t totalEmptyWork = 0;
    for (int i = 0; i < threadNum; i++) {
        totalEmptyWork += threadStats[i].emptyWork;
    }
    cout << "total empty work: " << totalEmptyWork << endl;

    // process cover
    if (!noverify) {
        vector<uintE> checkCover;
        checkCover.reserve(G.n);
        for (int i = 0; i < G.n; i++) {
            if (cover[i])
                checkCover.push_back(i);
        }
        if (!setcover::success<vertex>(G, checkCover)) abort();
    }
}

template<class vertex, bool usePrefetch>
void initialize(graph<vertex>& GA, commandLine P) {
    string algoType = P.getOptionValue("-type", "MQBucket");
    int threadNum = P.getOptionIntValue("-threads", 1);
    int queueNum = P.getOptionIntValue("-queues", 2);
    int bucketNum = P.getOptionLongValue("-buckets", 64);
    int batchSizePop = P.getOptionIntValue("-batch1", 1);
    int batchSizePush = P.getOptionIntValue("-batch2", 1);
    int stickiness = P.getOptionIntValue("-stick", 1);
    bool noverify = P.getOptionValue("-noverify");

    cout << "### Application: set-cover" << endl;
    cout << "### Graph: " << P.getArgument(0) << endl;
    cout << "### Workers: " << getWorkers() << endl;
    cout << "### n: " << GA.n << endl;
    cout << "### m: " << GA.m << endl;
    cout << "Type: " << algoType << endl;
    cout << "MQ thread num: " << threadNum << endl;
    cout << "MQ queue num: " << queueNum << endl;
    cout << "MQ batchsize pop: " << batchSizePop << endl;
    cout << "MQ batchsize push: " << batchSizePush << endl;
    cout << "MQ stickiness: " << stickiness << endl;
    cout << "MQ prefetch: " << usePrefetch << endl;

    // initialize to 0
    atomic<uint32_t>* cardinality = new atomic<uint32_t>[GA.n]();
    atomic<bool>* isElemCovered = new atomic<bool>[GA.n]();
    atomic<bool>* cover = new atomic<bool>[GA.n]();

    std::function<void(uint32_t)> prefetcher = [&] (uint32_t v) -> void {
        __builtin_prefetch(&cardinality[v], 0, 3);
    };

    if (algoType == "MQBucket") {
        cout << "MQ bucket num: " << bucketNum << endl;
        uintE m = 0;
        for (size_t s = 0; s < GA.n; s++) {
            m = max(m, GA.V[s].getOutDegree());
        }
        cout << "max cardinality = " << m << "\n";

        std::function<mbq::BucketID(uint32_t)> getBucketID = [&] (uint32_t v) -> mbq::BucketID {
            uint32_t card = cardinality[v].load(std::memory_order_seq_cst);
            return card;
        };
        using MQ_Bucket = mbq::MultiBucketQueue<
            decltype(getBucketID), decltype(prefetcher), 
            less<uintE>, uintE, uintE, usePrefetch>;
        MQ_Bucket wl(getBucketID, prefetcher, queueNum, threadNum, 0,
                 bucketNum, batchSizePop, batchSizePush, mbq::decreasing, stickiness, m);
        spawnTasks<vertex, MQ_Bucket>(GA, wl, threadNum, isElemCovered, cardinality, cover, noverify);

    } else if (algoType == "MQ") {
        using MQ = mbq::MultiQueue<decltype(prefetcher), less<PQElement>, uintE, uintE, usePrefetch>;
        MQ wl(prefetcher, queueNum, threadNum, batchSizePop, batchSizePush, stickiness);
        spawnTasks<vertex, MQ>(GA, wl, threadNum, isElemCovered, cardinality, cover, noverify);
    
    } else {
        cout << "Invalid type!\n";
    }
}

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
    bool usePrefetch = P.getOptionValue("-prefetch");
    if (usePrefetch) initialize<vertex, true>(GA, P);
    else initialize<vertex, false>(GA, P);
}
