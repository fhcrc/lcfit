.PHONY: all cmake-debug debug cmake-release release lcfit-compare lcfit-test test example clean doc

BUILD_DIR	:= _build
RELEASE_DIR	:= $(BUILD_DIR)/release
DEBUG_DIR 	:= $(BUILD_DIR)/debug

CMAKE_FLAGS	:= -D CMAKE_C_COMPILER=$(CC) -D CMAKE_CXX_COMPILER=$(CXX)

all: release

cmake-debug:
	mkdir -p $(DEBUG_DIR)
	cd $(DEBUG_DIR) && cmake $(CMAKE_FLAGS) -D CMAKE_BUILD_TYPE=Debug ../..

debug: cmake-debug
	$(MAKE) -C $(DEBUG_DIR)

cmake-release:
	mkdir -p $(RELEASE_DIR)
	cd $(RELEASE_DIR) && cmake $(CMAKE_FLAGS) -D CMAKE_BUILD_TYPE=Release ../..

release: cmake-release
	$(MAKE) -C $(RELEASE_DIR)

lcfit-compare: debug
	$(MAKE) -C $(DEBUG_DIR) $@

lcfit-test: debug
	$(MAKE) -C $(DEBUG_DIR) $@

test: lcfit-test
	$(DEBUG_DIR)/test/lcfit-test 2> lcfit-test.log

example: lcfit-compare
	$(MAKE) -C example

clean:
	rm -rf $(BUILD_DIR)

doc:
	doxygen
