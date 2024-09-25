ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

ifdef PD
PD = -DPD
endif

ifdef BYTE
CODE = -DBYTE
else ifdef NIBBLE
CODE = -DNIBBLE
else
CODE = -DBYTERLE
endif

# Use clang to compile
PCC = clang++
PCFLAGS = -std=c++17 -fopencilk -O3 -DCILK $(INTT) $(INTE) $(CODE) $(PD) $(FLAGS)
PLFLAGS = 

COMMON= ligra.h edgeMap_utils.h index_map.h sequence.h maybe.h binary_search.h graph.h compressedVertex.h vertex.h utils.h IO.h parallel.h gettime.h quickSort.h blockRadixSort.h transpose.h parseCommandLine.h byte.h byteRLE.h nibble.h byte-pd.h byteRLE-pd.h nibble-pd.h vertexSubset.h encoder.C bucket.h counting_sort.h dyn_arr.h edgeMapReduce.h histogram.h sequentialHT.h

CPS= MultiQueue.h MultiBucketQueue.h BucketStructs.h

BKTMQ= SetCover_MQ

JUL= SetCover DeltaStepping

ALL= DeltaStepping SetCover

SetCover: SetCover.C $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

DeltaStepping: DeltaStepping.C $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

SetCover_MQ: SetCover_MQ.C $(COMMON) $(CPS)
	clang++ -std=c++17 -O3 -pthread $(INTT) $(INTE) $(CODE) $(PD) -o $@ $<

all: $(BKTMQ) $(JUL)

$(COMMON):
	ln -s ../../ligra/$@ .

$(CPS):
	ln -s ../../../CPS/include/$@ .

.PHONY : clean

clean :
	rm -f *.o $(BKTMQ) $(JUL)

cleansrc :
	rm -f *.o $(BKTMQ) $(JUL)
	rm $(COMMON)
	rm $(CPS)

