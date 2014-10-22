# gcc for c compilation and g++ for cpp compilation
CC       = gcc
CXX      = g++



# include some libraries
INCLUDE  = 



# compilation options: -g for debugging info and -Wall for enable all warnings
CFLAGS   = -g -Wall $(INCLUDE)
CXXFLAGS = -g -Wall $(INCLUDE)



# linking options: -g for debugging info
LDFLAGS  = -g



# list the libraries that we need to link with in LDLIBS
# eg: use "-lm" for the math library, if we use sqrt() function.
LDLIBS   = -lm



# the 1st target gets built when we type "make", which is our executable.
# when we don't specify the linking rule, instead, we rely on one of make's impilict rules:
#
# $(CC)$(LDFLAGS) <all-dependent-.o-file> $(LDLIBS)
#
# alsonot that make assumes that main depends on main.o which means we could omit it.

main: CharacteristicFuncLib.o FourierLib.o OptionPricing.o main.o

# because main depends on those CharacteristicFunclib.h, FourierLib.h, OptionPricing.h
# main.o will get recompiled if either of the header files get modified
# make already know main.o depends on main.c, so we don't have to say it. make use
# following impilicit rule to compile a .c file into a .o file:
# 
# $(CC) -c $(CFLAGS) <the-.c-file>
#

main.o: CharacteristicFuncLib.h FourierLib.h OptionPricing.h main.cpp


# and Charac.o, Fourier.o, Option.o depend on their own .h files, again, make knows then also depend on the .cpp files.

CharacteristicFuncLib.o:  CharacteristicFuncLib.h CharacteristicFuncLib.cpp
FourierLib.o           :  FourierLib.h FourierLib.cpp
OptionPricing.o        :  OptionPricing.h OptionPricing.cpp



# always provide the "clean" target that removes intermediate files
# what we remove depend our own choice of coding tools, different editors generate
# different backup files for example
# "clean" target is not a file name, we sold tell make it's a "phony" target.

.PHONY: clean
clean:
	rm -f *.o a.out core main *~


# "all" target is useful if my makefile builds multiple programs.
# usually we make it first do "clean", and re build the main target.
.PHONY: all
all:
	clean main

  

