.PHONY: all lcfit-compare lcfit-test setup-cmake clean run test doc release debug

BUILD := _build

all: release

release: CMAKE_BUILD_TYPE=Release
release: BUILD_DIR=$(BUILD)/release
release: lcfit-compare test

debug: CMAKE_BUILD_TYPE=Debug
debug: BUILD_DIR=$(BUILD)/debug
debug: lcfit-compare test

data.pdf: data.csv
	Rscript plot_fits.R data.csv data.maxima.csv data.fit.csv data.pdf

data.csv: lcfit-compare
	$(BUILD_DIR)/lcfit-compare \
		param=data/test.params.bpp

lcfit-compare: setup-cmake
	+make -C$(BUILD_DIR) $@

lcfit-test: setup-cmake
	+make -C$(BUILD_DIR) $@

test: lcfit-test
	$(BUILD_DIR)/$<

setup-cmake:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) ../..

clean:
	rm -rf $(BUILD)

doc:
	doxygen

style:
	astyle  -A3 \
	        --pad-oper \
	        --unpad-paren \
	        --keep-one-line-blocks \
	        --keep-one-line-statements \
	        --suffix=none \
	        --formatted \
	        --lineend=linux \
					--align-pointer=type \
	        `find src -regextype posix-extended -regex ".*\.(cc|h|hpp)$$"`
