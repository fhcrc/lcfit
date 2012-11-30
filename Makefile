.PHONY: all lcfit-compare lcfit-test setup-cmake clean run

BUILD := _build

run: lcfit-compare
	$(BUILD)/lcfit-compare data/test.tre data/test.fasta
	graph -T svg < data.dat > data.svg

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
