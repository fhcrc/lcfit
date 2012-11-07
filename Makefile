LIBS =

CC ?= gcc
CXX ?= g++

CFLAGS := -g -I. $(CFLAGS)
LFLAGS := $(LFLAGS) -lm -lgsl -lgslcblas

all: run

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

test: test.o quintic_C.o
	$(CC) -o $@ $< quintic_C.o -g $(LFLAGS) $(CFLAGS)

run: test
	./test

style:
	astyle  -A3 \
		--pad-oper \
		--unpad-paren \
		--keep-one-line-blocks \
		--keep-one-line-statements \
		--suffix=none \
		--formatted \
		--lineend=linux \
		`find . -regextype posix-extended -regex ".*\.(c|h)$$"`
