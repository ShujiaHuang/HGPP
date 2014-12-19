BAMTOOLS_ROOT = ../include/bamtools
CXXFLAGS+=-Wl,-rpath,$(BAMTOOLS_ROOT)/lib -Wall -O3 -I $(BAMTOOLS_ROOT)/include -L $(BAMTOOLS_ROOT)/lib
CXXFLAGS+=-I../include/utils

bin = ../../bin
obj = utility.o sam.o
anchor_fosmid: $(obj) anchor_fosmid.cpp anchor_fosmid.h
	$(CXX) $(CXXFLAGS) anchor_fosmid.cpp anchor_fosmid.h $(obj) -lbamtools -o $(bin)/anchor/anchor_fosmid -lz

utility.o: ../include/utils/utility.cpp ../include/utils/utility.h
	$(CXX) -c -o utility.o ../include/utils/utility.cpp

sam.o: ../include/utils/sam.cpp ../include/utils/sam.h
	$(CXX) -c -o sam.o ../include/utils/sam.cpp

clean:
	rm -f $(obj)

.PHONY: clean

