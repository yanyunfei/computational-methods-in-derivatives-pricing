CC  = clang++
CXX = clang++

INCLUDES = 
CFLAGS   = -g -Wall $(INCLUDES)
CXXflags = -g -Wall $(INCLUDES)

LDFLAGS  = -g
LDLIBS   = -g

heston:
	g++ heston.cpp -lgsl -lgslcblas -o heston	


.PHONY: clean
clean:
	rm -rf *.c~ *.txt~ a.out *.cpp~ Makefile~ heston 1 *.cc~ *.txt~
.PHONY: all
all:
	clean a.out
