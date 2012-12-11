#!/bin/sh

#nestagg delim -d runs lcfit_plots.txt -t -o lcfit_plots_agg.txt -k initial,ml_tree,model_name
nestagg delim -d runs lcfit_maxima.csv -o lcfit_maxima.csv -k ml_tree,model_name
./aggregate.R
