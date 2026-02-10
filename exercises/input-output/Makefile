CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
# LINK.o = $(CXX) $(LDFLAGS)
LDLIBS += -lstdc++ -lm

all: cmdline stdin filestream


cmdline : cmdline.o
	g++ cmdline.o -o cmdline
cmdline.o : cmdline.cc
	g++ -c cmdline.cc -o cmdline.o


stdin : stdin.o
	g++ stdin.o -o stdin
stdin.o : stdin.cc
	g++ -c stdin.cc -o stdin.o


filestream : filestream.o
	g++ filestream.o -o filestream
filestream.o : filestream.cc
	g++ -c filestream.cc -o filestream.o


.PHONY:clean
clean:
	$(RM) cmdline stdin filestream *.o

test:
	@echo $(CXX)
	@echo $(RM)
