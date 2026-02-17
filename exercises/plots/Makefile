CXXFLAGS += -Wall -std=c++23
LDLIBS += -lstdc++ -lm

fig1.svg: fig1.gpi data.data
	gnuplot $<

data.data: main Makefile
	./$< > $@

.PHONY:clean
clean:
	$(RM) main