.PHONY: all lcfit-compare lcfit-test setup-cmake clean run

BUILD := _build


run: lcfit-compare data.pdf data.csv

data.pdf: plot_fits.R data.csv
	Rscript plot_fits.R data.csv data.pdf

data.csv: lcfit-compare
	$(BUILD)/lcfit-compare data/test.tre data/test.fasta

all: lcfit-compare

lcfit-compare: setup-cmake
	+make -C$(BUILD) $@

setup-cmake:
	mkdir -p $(BUILD)
	cd $(BUILD) && cmake -DCMAKE_BUILD_TYPE=Debug ..

clean:
	rm -rf $(BUILD)

style:
	astyle  -A3 \
	        --pad-oper \
	        --unpad-paren \
	        --keep-one-line-blocks \
	        --keep-one-line-statements \
	        --suffix=none \
	        --formatted \
	        --lineend=linux \
	        `find src -regextype posix-extended -regex ".*\.(cc|h|hpp)$$"`
