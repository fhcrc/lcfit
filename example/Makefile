likelihoods.pdf: likelihoods.csv
	Rscript plot_fits.R likelihoods.csv likelihoods.maxima.csv likelihoods.fit.csv $@

likelihoods.csv:
	../_build/release/lcfit_cpp_src/lcfit-compare \
		param=test.params.bpp
