CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt timing.svg

Out.txt timing.data : main ; ./$< > $@
main.o : matrix.h qr.h
timing.svg : timing.gpi timing.data
	gnuplot $<

clean:
	$(RM) main *.o *.txt *.data *.svg
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)