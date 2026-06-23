CXX = g++
CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt timing.data wavefunctions.data timing.svg wavefunctions.svg

Out.txt: main 
	./main > Out.txt
	
# main: main.o 
# main.o: main.cc root.h run_timing.h

main: main.o run_timing.o hamiltonian.o
main.o: main.cc root.h
run_timing.o: run_timing.cc run_timing.h
hamiltonian.o: hamiltonian.cc hamiltonian.h

timing.svg : timing.gpi timing.data
	gnuplot $<
wavefunctions.svg : wavefunctions.gpi wavefunctions.data
	gnuplot $<

clean:
	$(RM) main *.o *.txt *.svg
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)