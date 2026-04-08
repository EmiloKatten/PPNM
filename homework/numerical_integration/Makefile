CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt erf.svg erf_acc.svg

Out.txt data_erf.txt data_erf_acc.txt: main 
	./main > Out.txt
	
main : main.o integrator.o
main.o : main.cc
integrator.o : integrator.cc integrator.h

erf.svg : erf.gpi data_erf.txt
	gnuplot $<

erf_acc.svg : erf_acc.gpi data_erf_acc.txt
	gnuplot $<

clean:
	$(RM) main *.o *.txt *.svg
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)