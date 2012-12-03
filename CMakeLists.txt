cmake_minimum_required (VERSION 2.8.1)

project(lcfit CXX)

set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -pedantic")

# STS
#add_subdirectory(sts)
#include_directories(sts/lib/include)

set(LCFIT_COMPARE_CPP ${CMAKE_CURRENT_SOURCE_DIR}/src/compare.cc)
#set(LCFIT_TEST_CPP ${CMAKE_CURRENT_SOURCE_DIR}/src/test.cc)
include_directories(src)

set(LCFIT_LIBS
  bpp-core
  bpp-seq
  bpp-phyl
  gsl
  gslcblas)

# LCFIT
add_executable(lcfit-compare ${LCFIT_COMPARE_CPP})
target_link_libraries(lcfit-compare ${LCFIT_LIBS})
