CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt ann.svg

Out.txt data.txt fit.txt integral.txt: main 
	./main > Out.txt
	
main : main.o 
main.o : main.cc ann.h

ann.svg : ann.gpi data.txt fit.txt integral.txt
	gnuplot $<

clean:
	$(RM) main *.o *.txt *.svg log
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)