CXXFLAGS+=-I../include/utils -Wall

bin = ../../bin
obj = ErrCorr.o utility.o
ErrCorr: $(obj)
	$(CXX) $(LDFLAGS) -o $(bin)/phasing/ErrCorr $(obj)

utility.o: ../include/utils/utility.cpp
	$(CXX) -c -o utility.o ../include/utils/utility.cpp

ErrCorr.o: ErrCorr.cpp
	$(CXX) $(CXXFLAGS) -c -o ErrCorr.o ErrCorr.cpp


clean:
	rm -f $(obj)

.PHONY: clean
