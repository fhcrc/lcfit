LIBS =

CC ?= gcc
CXX ?= g++

CFLAGS := -g -I. $(CFLAGS) -std=c++0x
LFLAGS := $(LFLAGS) -lm -lgsl -lgslcblas -lbpp-core -lbpp-seq -lbpp-phyl

all: run

_build:
	mkdir -p _build

_build/%.o: src/%.cc src/lcfit.h _build
	$(CXX) -c -o $@ $< $(CFLAGS)

compare: _build/compare.o
	$(CXX) -o $@ $< -g $(LFLAGS) $(CFLAGS)

test: _build/test.o
	$(CXX) -o $@ $< -g $(LFLAGS) $(CFLAGS)

run: compare
	./compare
	graph -T svg < data.dat > data.svg

clean:
	rm -rf _build compare test

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
