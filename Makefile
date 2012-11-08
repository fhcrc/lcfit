LIBS =

CC ?= gcc
CXX ?= g++

CFLAGS := -g -I. $(CFLAGS)
LFLAGS := $(LFLAGS) -lm -lgsl -lgslcblas

all: run

%.o: %.c lcfit.h
	$(CC) -c -o $@ $< $(CFLAGS)

test: test.o
	$(CC) -o $@ $< -g $(LFLAGS) $(CFLAGS)

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
