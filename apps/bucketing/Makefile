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

OPT = -O3 -mtune=native

#compilers
ifdef CILK
PCC = clang++
PCFLAGS = -std=c++14 -fopencilk $(OPT) -DOPENCILK $(INTT) $(INTE) $(CODE) $(PD) $(FLAGS)
PLFLAGS = -fopencilk

else ifdef MKLROOT
PCC = icpc
PCFLAGS = -std=c++14 $(OPT) -DCILKP $(INTT) $(INTE) $(CODE) $(PD)

else ifdef OPENMP
PCC = g++
PCFLAGS = -std=c++14 -fopenmp -march=native $(OPT) -DOPENMP $(INTT) $(INTE) $(CODE) $(PD)

else
PCC = g++
PCFLAGS = -std=c++17 $(OPT)
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
