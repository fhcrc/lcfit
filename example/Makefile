.PHONY: all

all: curves.pdf

lnl.csv:
	../_build/debug/lcfit_cpp_src/lcfit-compare param=test.params.bpp

curves.pdf: lnl.csv
	Rscript plot_example.R
