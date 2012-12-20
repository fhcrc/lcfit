cmake_minimum_required (VERSION 2.8.1)

# Build the LCFIT C library
project(libclfit C)
set(CMAKE_C_FLAGS "-std=c99 -Wall -pedantic")
set(LCFIT_LIB_C_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/lcfit.c)
add_library(lcfit-static STATIC ${LCFIT_LIB_C_FILES})

set(LCFIT_TEST_LIBS
  lcfit-static
  check
  gsl
  gslcblas)

# Test
include_directories(src)
set(LCFIT_TEST_C_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/test/check_basic.c)
add_executable(lcfit-test ${LCFIT_TEST_C_FILES})
target_link_libraries(lcfit-test ${LCFIT_TEST_LIBS})


# Comparison
project(lcfit CXX)
set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -pedantic")

set(LCFIT_COMPARE_CPP ${CMAKE_CURRENT_SOURCE_DIR}/src/compare.cc)
include_directories(src)

set(LCFIT_LIBS
  bpp-core
  bpp-seq
  bpp-phyl
  gsl
  gslcblas)

# LCFIT
add_executable(lcfit-compare ${LCFIT_COMPARE_CPP})
target_link_libraries(lcfit-compare lcfit-static ${LCFIT_LIBS})


