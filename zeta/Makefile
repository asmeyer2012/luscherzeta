SRCS = $(wildcard *.cc)
PROGS = $(patsubst main%.cc,exe%,$(SRCS))
CC = g++
## get from output of $(gsl-config --cflags --libs) $(pkg-config --cflags --libs blas)
LIBS = -lgsl
INC =

## compiles with g++ version 6.1.0, 6.3.0
## does not compile with g++ version 4.8.5

all: libczeta.so $(PROGS)

exe: main.cc
	${CC} -Wall -O3 -std=c++11 -o $@ $< ${INC} ${LIBS}

libczeta.so: czeta.cc
	${CC} -Wall -shared -O3 -fPIC -std=c++11 -o $@ $< ${INC} ${LIBS}

clean:
	rm -rf exe* *.o *.so

