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

#compilers
ifdef CILK
PCC = /nfs/ug/groups/yangste9/cilk/build/bin/clang++
PCFLAGS = -std=c++14 -fopencilk -O3 -DCILK $(INTT) $(INTE) $(CODE) $(PD) $(FLAGS)
PLFLAGS = -fopencilk

else ifdef MKLROOT
PCC = icpc
PCFLAGS = -std=c++14 -O3 -DCILKP $(INTT) $(INTE) $(CODE) $(PD)

else ifdef OPENMP
PCC = g++
PCFLAGS = -std=c++14 -fopenmp -march=native -O3 -DOPENMP $(INTT) $(INTE) $(CODE) $(PD)

else
PCC = /nfs/ug/homes-0/y/yangste9/gcc11/bin/g++
PCFLAGS = -std=c++17 -O3
endif

COMMON= ligra.h edgeMap_utils.h index_map.h sequence.h maybe.h binary_search.h graph.h compressedVertex.h vertex.h utils.h IO.h parallel.h gettime.h quickSort.h blockRadixSort.h transpose.h parseCommandLine.h byte.h byteRLE.h nibble.h byte-pd.h byteRLE-pd.h nibble-pd.h vertexSubset.h encoder.C bucket.h counting_sort.h dyn_arr.h edgeMapReduce.h histogram.h sequentialHT.h

ALL= DeltaStepping KCore KCoreSerial SetCover

all: $(ALL)

% : %.C $(COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

$(COMMON):
	ln -s ../../ligra/$@ .

.PHONY : clean

clean :
	rm -f *.o $(ALL)

cleansrc :
	rm -f *.o $(ALL)
	rm $(COMMON)
