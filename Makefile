LIBS =

CC ?= gcc
CXX ?= g++

CFLAGS := -g -I. $(CFLAGS)
LFLAGS := $(LFLAGS) -lm

all: run

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

test: test.o quintic_C.o
	$(CC) -o $@ $< quintic_C.o -g $(LFLAGS) $(CFLAGS)

run: test
	./test

clean:
	rm -f *.o test
