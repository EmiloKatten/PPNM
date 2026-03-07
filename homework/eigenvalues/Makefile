CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

Out.txt:main ; ./$< -rmax 10 -dr 0.3 > $@
main.o: matrix.h evd.h
dr.svg: dr.gpi data_dr.txt
	gnuplot $<
rmax.svg: rmax.gpi data_rmax.txt
	gnuplot $<
eigenfunctions.svg: eigenfunctions.gpi eigenfunctions.txt
	gnuplot $<
clean:
	$(RM) main *.o
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)