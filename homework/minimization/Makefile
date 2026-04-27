CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt higgs.svg

Out.txt fit.txt: main 
	./main < higgs.dat 1> Out.txt
	
main : main.o 
main.o : main.cc newton.h

higgs.svg : higgs.gpi higgs.dat fit.txt
	gnuplot $<

clean:
	$(RM) main *.o *.txt *.svg log
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)