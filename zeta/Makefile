SRCS = $(wildcard *.cc)
PROGS = $(patsubst main%.cc,main.o%,$(SRCS))
CC = g++
LIBS = -lgsl -L/usr/lib/libblas -lblas -lgslcblas
INC =

## compiles with g++ version 6.1.0, 6.3.0
## does not compile with g++ version 4.8.5

all: libczeta.so $(PROGS)

main.o: main.cc
	${CC} -Wall -O3 -std=c++11 ${INC} -o $@ $^ ${LIBS}

libczeta.so: czeta.cc
	${CC} -Wall -shared -O3 -fPIC -std=c++11 ${INC} -o $@ $^ ${LIBS}

clean:
	rm -rf exe *.o *.so

