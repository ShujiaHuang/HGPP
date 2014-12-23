BAMTOOLS_ROOT = $(realpath ../../src/include/bamtools)
CXXFLAGS+=-Wl,-rpath,$(BAMTOOLS_ROOT)/lib -Wall -O3 -I $(BAMTOOLS_ROOT)/include -L $(BAMTOOLS_ROOT)/lib
CXXFLAGS+=-I../../src/include/utils

bin = ./
obj = utility.o sam.o
DistinguishReads: $(obj) DistinguishReads.cpp anchor_variation.h
	$(CXX) $(CXXFLAGS) DistinguishReads.cpp anchor_variation.h $(obj) -lbamtools -o $(bin)/DistinguishReads -lz

utility.o: ../../src/include/utils/utility.cpp
	$(CXX) -c -o utility.o ../../src/include/utils/utility.cpp

sam.o: ../../src/include/utils/sam.cpp ../../src/include/utils/sam.h 
	$(CXX) -c -o sam.o ../../src/include/utils/sam.cpp

clean:
	rm -f $(obj)

.PHONY: clean

