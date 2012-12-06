#!/bin/sh

nestagg delim -d runs lcfit_plots.txt -t -o lcfit_plots_agg.txt -k initial,ml_tree
./aggregate.R
