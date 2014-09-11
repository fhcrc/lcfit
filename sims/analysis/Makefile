
%.html : %.Rmd
	R --slave -e 'rmarkdown::render("$(firstword $^)")'

fitting.html : fitting.Rmd

fitting.html : utils.R

aggfit.html : aggfit.Rmd

aggfit.html : aggutils.R

# Cribbed from http://dirk.eddelbuettel.com/code/rcpp.html  
# "Usage for package building"
#
PKG_CXXFLAGS := -I$(shell Rscript -e "cat(R.home('include'))") $(shell Rscript -e "Rcpp:::CxxFlags()")
PKG_LIBS := $(shell Rscript -e "cat(R.home('lib'))") $(shell Rscript -e "Rcpp:::LdFlags()")

CPPFLAGS =  $(PKG_CXXFLAGS) -fpic
CFLAGS =  $(PKG_CXXFLAGS) -fpic
LDLIBS = $(patsubst %,-L%,$(PKG_LIBS))


lcfit_rcpp.o : lcfit_rcpp.cpp


lcfit.o : ../lcfit_src/lcfit.c 
	$(COMPILE.c) $(OUTPUT_OPTION) $<

lcfit_rcpp.so : lcfit_rcpp.o lcfit.o
	g++ -shared -L/usr/local/lib -o $@ lcfit_rcpp.o lcfit.o $(LDLIBS) -lR -lgsl

clean:
	rm -f lcfit_rcpp.o lcfit.o
