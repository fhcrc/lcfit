# To build lcfit-compare, you need a specific version of bpp in your path
#
#	module load bpp/master-201404114
#
# To build a debug version of lcfit,
#
#	make -C _build/debug/ lcfit-compare
#
# This will leave lcfit-compare in ./_build/debug/lcfit_cpp_src/lcfit-compare



.PHONY: all lcfit-compare lcfit-test setup-cmake clean run test doc release debug build-all example

BUILD := _build
CMAKE_BUILD_TYPE = Release
BUILD_DIR = $(BUILD)/release

all: release

release: CMAKE_BUILD_TYPE=Release
release: BUILD_DIR=$(BUILD)/release
release: build-all lcfit-r

debug: CMAKE_BUILD_TYPE=Debug
debug: BUILD_DIR=$(BUILD)/debug
debug: build-all

build-all: setup-cmake
	$(MAKE) -C $(BUILD_DIR)

test: lcfit-test

example: release lcfit-compare
example:
	$(MAKE) -C example

lcfit-r:
	R -e "Rcpp::compileAttributes('lcfit_R')"
	R CMD INSTALL --library=$(PWD)/sims/venv/lib/R lcfit_R

# lcfit-compare: release
lcfit-compare: setup-cmake
	$(MAKE) -C$(BUILD_DIR) $@

lcfit-test: CMAKE_BUILD_TYPE=Debug
lcfit-test: BUILD_DIR=$(BUILD)/debug
lcfit-test: setup-cmake
	$(MAKE) -C$(BUILD_DIR) $@
	$(BUILD_DIR)/lcfit_cpp_src/lcfit-test -s

setup-cmake:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) ../..

clean:
	rm -rf $(BUILD)
	R CMD REMOVE --library=$(PWD)/sims/venv/lib/R lcfit

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
