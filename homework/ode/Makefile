CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt oscillator.svg planets.svg

Out.txt data_oscillator.txt data_planets1.txt data_planets2.txt data_planets3.txt: main 
	./main > Out.txt
	
main : main.o
main.o : main.cc matrix.h rkstep.h

oscillator.svg : oscillator.gpi data_oscillator.txt
	gnuplot $<
planets.svg : planets.gpi data_planets1.txt data_planets2.txt data_planets3.txt
	gnuplot $<

clean:
	$(RM) main *.o *.txt *.svg
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)