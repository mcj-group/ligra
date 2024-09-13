#!/bin/bash

input_dir=../../inputs_pbbs/ag

graphs=(USA LJ TW HYPERH)

stickinessbucket=(8 4 8 1)
stickinessheap=(8 4 8 1)
batches1=(256 256 32 256)
batches2=(256 256 256 256)

# MQBucket
for i in {0..3}
do
	echo "Results" > results_sc_48_${graphs[i]}_MQBucket.log
	./SetCover_MQ -rounds 10 -delta 0 -threads 48 -queues 192 -buckets 64 -stick ${stickinessbucket[i]} -batch1 256 -batch2 256 -type MQBucket -prefetch ${input_dir}/${graphs[i]}.ag  >> results_sc_48_${graphs[i]}_MQBucket.log
done

for i in {0..3}
do
        echo "Results" > results_sc_1_${graphs[i]}_MQBucket.log
        ./SetCover_MQ -rounds 10 -delta 0 -threads 1 -queues 4 -buckets 64 -stick 1 -batch1 256 -batch2 256 -type MQBucket -prefetch ${input_dir}/${graphs[i]}.ag  >> results_sc_1_${graphs[i]}_MQBucket.log
done


# MQ
for i in {0..3}
do
        echo "Results" > results_sc_48_${graphs[i]}_MQ.log
        ./SetCover_MQ -rounds 10 -delta 0 -threads 48 -queues 192 -buckets 64 -stick ${stickinessheap[i]} -batch1 ${batches1[i]} -batch2 ${batches2[i]} -type MQ -prefetch ${input_dir}/${graphs[i]}.ag  >> results_sc_48_${graphs[i]}_MQ.log
done

for i in {0..3}
do
        echo "Results" > results_sc_1_${graphs[i]}_MQ.log
        ./SetCover_MQ -rounds 10 -delta 0 -threads 1 -queues 4 -buckets 64 -stick 1 -batch1 ${batches1[i]} -batch2 ${batches2[i]} -type MQ -prefetch ${input_dir}/${graphs[i]}.ag  >> results_sc_1_${graphs[i]}_MQ.log
done

# MQPlain
for i in {0..3}
do
    echo "Results" > results_sc_48_${graphs[i]}_MQPlain.log
    ./SetCover_MQ -rounds 10 -threads 48 -queues 192 -batch1 1 -batch2 1 -type MQ ${input_dir}/${graphs[i]}.ag  >> results_sc_48_${graphs[i]}_MQPlain.log
done

for i in {0..3}
do
    echo "Results" > results_sc_1_${graphs[i]}_MQPlain.log
    ./SetCover_MQ -rounds 10 -threads 1 -queues 4 -batch1 1 -batch2 1 -type MQ ${input_dir}/${graphs[i]}.ag  >> results_sc_1_${graphs[i]}_MQPlain.log
done

