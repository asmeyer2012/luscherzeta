SRCS = $(wildcard *.cc)
PROGS = $(patsubst main%.cc,exe%,$(SRCS))
CC = g++
LIBS = -lgsl -L/usr/lib/libblas -lblas
INC =

all: libczeta.so $(PROGS)

exe: main.cc
	${CC} -Wall -O3 -std=c++11 ${INC} -o $@ ${LIBS} $<

libczeta.so: czeta.cc
	${CC} -Wall -shared -O3 -fPIC -std=c++11 ${INC} -o $@ ${LIBS} $<

clean:
	rm -rf exe *.o *.so

