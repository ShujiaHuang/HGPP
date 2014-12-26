CXXFLAGS+=-I../include/utils -Wall

bin = ../../bin
obj = cnn.o utility.o
cnn: $(obj)
	$(CXX) $(LDFLAGS) -o $(bin)/phasing/cnn $(obj)

utility.o: ../include/utils/utility.cpp
	$(CXX) -c -o utility.o ../include/utils/utility.cpp

cnn.o: cnn.cpp
	$(CXX) $(CXXFLAGS) -c -o cnn.o cnn.cpp


clean:
	rm -f $(obj)

.PHONY: clean
