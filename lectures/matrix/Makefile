CXXFLAGS = -std=c++23 -Wall -O
LDLIBS = -lstdc++ -lm

Out.txt: main
	./$< > $@

main.o: matrix.h

clean:
	$(RM) *.o main