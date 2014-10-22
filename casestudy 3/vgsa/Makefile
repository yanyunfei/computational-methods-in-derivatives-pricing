CC  = gcc
CXX = g++

INCLUDES = 
CFLAGS   = -g -Wall $(INCLUDES)
CXXflags = -g -Wall $(INCLUDES)

LDFLAGS  = -g
LDLIBS   = -g

vgsa:
	g++ vgsa.cc -lgsl -lgslcblas -o vgsa	


.PHONY: clean
clean:
	rm -rf *.c~ *.txt~ a.out *.cpp~ Makefile~ heston 1 *.cc~ vgsa
.PHONY: all
all:
	clean a.out
