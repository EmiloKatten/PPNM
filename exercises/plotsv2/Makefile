CXXFLAGS += -Wall -std=c++23
LDLIBS += -lstdc++ -lm

all: erf.svg sgamma.svg lngamma.svg

erf.svg: erf.gpi erf.data
	gnuplot $<
sgamma.svg: sgamma.gpi sgamma.data
	gnuplot $<
lngamma.svg: lngamma.gpi lngamma.data
	gnuplot $<

erf.data: main Makefile
	./$< erf > $@
sgamma.data: main Makefile
	./$< sgamma > $@
lngamma.data: main Makefile
	./$< lngamma > $@

.PHONY:clean
clean:
	$(RM) main