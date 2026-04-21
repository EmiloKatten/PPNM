CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt wavefunction.svg

Out.txt wavefunction.txt: main 
	./main > Out.txt
	
main : main.o 
main.o : main.cc root.h rkstep.h

wavefunction.svg : wavefunction.gpi wavefunction.txt
	gnuplot $<

clean:
	$(RM) main *.o *.txt *.svg
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)