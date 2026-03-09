CXXFLAGS += -std=c++23 -Wall -Wextra #-DLONG_DOUBLE
LDLIBS += -lstdc++ -lm

default: Out.txt

Out.txt : main ; ./$< > $@
main.o : matrix.h qr.h
clean:
	$(RM) main *.o *.txt
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)