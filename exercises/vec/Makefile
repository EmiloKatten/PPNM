CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
# LINK.o = $(CXX) $(LDFLAGS)
LDLIBS += -lstdc++ -lm

all: Out.txt
	cat $<

Out.txt : main
	./main > Out.txt

main : main.o vec.o 
main.o : vec.h main.cc
vec.o: vec.h vec.cc


.PHONY:clean
clean:
	$(RM) main main.o Out.txt *.o

test:
	@echo $(CXX)
	@echo $(RM)
