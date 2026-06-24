CXX = g++
CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

.PHONY: default clean

default: Out.txt timing.svg comp_timing.svg wavefunctions.svg perturbation.svg

Out.txt timing.data comp_timing.data wavefunctions.data perturbation.data: main 
	./main > Out.txt

main: main.o run_timing.o hamiltonian.o

main.o: main.cc root.h
run_timing.o : run_timing.cc run_timing.h
hamiltonian.o : hamiltonian.cc hamiltonian.h

%.svg: %.gpi %.data
	gnuplot $<

clean:
	$(RM) main *.o *.txt *.data *.svg