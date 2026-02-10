CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
# LINK.o = $(CXX) $(LDFLAGS)
LDLIBS += -lstdc++ -lm

all: Out.txt
	cat $<

Out.txt : main
	./main > Out.txt

main : main.o machine.o addition.o compare.o
main.o : machine.h addition.h compare.h main.cc
machine.o: machine.h machine.cc
addition.o: addition.h addition.cc
compare.o: compare.h compare.cc


.PHONY:clean
clean:
	$(RM) main main.o *.o

test:
	@echo $(CXX)
	@echo $(RM)
