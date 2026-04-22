CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt error.svg

Out.txt data_error.txt : main 
	./main > Out.txt
	
main : main.o 
main.o : main.cc mc.h lcg.h

error.svg : error.gpi data_error.txt quasi_error.txt
	gnuplot $<


clean:
	$(RM) main *.o *.txt *.svg
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)