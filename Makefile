SRCS = $(wildcard *.cc)
PROGS = $(patsubst main%.cc,exe%,$(SRCS))
CC = g++
#LIBS =
#LIBS = -lgmp -lgsl
#LIBS = -lgsl
#LIBS = -lgsl -lblas
#LIBS = -lgsl -L/usr/lib -lblas
#LIBS = -lgsl -L/usr/lib/ -lblas
#LIBS = -lgsl -L/usr/lib/libblas -lblas
LIBS = -lgsl -L/usr/lib/libblas -lblas
INC =

all: $(PROGS)

exe: main.cc
	${CC} -O3 -std=c++11 ${INC} -o $@ ${LIBS} $<

clean:
	rm -rf exe *.o

#
#-mcmodel=large

