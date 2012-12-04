#!/bin/sh

nestagg delim -d runs lcfit_plots.txt -t -o lcfit_plots_agg.txt -k initial
./aggregate.R
