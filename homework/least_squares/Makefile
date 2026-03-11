CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt data.svg 

Out.txt data.txt fit.txt: main ; ./$< > $@
main.o : matrix.h qr.h ls.h
data.svg : data.gpi data.txt fit.txt
	gnuplot $<
clean:
	$(RM) main *.o *.txt *.svg
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)