CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: linterp.svg qinterp_f.svg qinterp_g.svg

ldata.txt linterp.txt qdata.txt qinterp.txt : main 
	./main
main : main.o qspline.o
main.o : main.cc qspline.h
qspline.o : qspline.cc qspline.h

linterp.svg : linterp.gpi ldata.txt linterp.txt
	gnuplot $<
qinterp_f.svg : qinterp_f.gpi qdata.txt qinterp.txt
	gnuplot $<
qinterp_g.svg : qinterp_g.gpi qdata.txt qinterp.txt
	gnuplot $<

clean:
	$(RM) main *.o *.txt *.svg
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)