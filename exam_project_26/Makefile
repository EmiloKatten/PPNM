CXX = g++
CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt  timing.svg comp_timing.svg wavefunctions.svg perturbation.svg

Out.txt timing.data wavefunctions.data test.data: main 
	./main > Out.txt
	
# main: main.o 
# main.o: main.cc root.h run_timing.h

main: main.o run_timing.o hamiltonian.o
main.o: main.cc root.h
run_timing.o: run_timing.cc run_timing.h
hamiltonian.o: hamiltonian.cc hamiltonian.h

timing.svg : timing.gpi timing.data
	gnuplot $<
comp_timing.svg : comp_timing.gpi comp_timing.data
	gnuplot $<
wavefunctions.svg : wavefunctions.gpi wavefunctions.data
	gnuplot $<
perturbation.svg : perturbation.gpi perturbation.data
	gnuplot $<
	
clean:
	$(RM) main *.o *.txt *.data *.svg
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)