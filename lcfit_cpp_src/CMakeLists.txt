cmake_minimum_required (VERSION 2.8.1)

# Comparison
project(lcfitcompare C CXX)
set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -pedantic")

set(LCFIT_COMPARE_CPP
  ${CMAKE_CURRENT_SOURCE_DIR}/lcfit_compare.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lcfit_cpp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../lcfit_src/lcfit.c)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lcfit_src)

set(LCFIT_COMPARE_LIBS
  bpp-core
  bpp-seq
  bpp-phyl
  gsl
  gslcblas
  lcfit-static)

# LCFIT
add_executable(lcfit-compare EXCLUDE_FROM_ALL ${LCFIT_COMPARE_CPP})
target_link_libraries(lcfit-compare ${LCFIT_COMPARE_LIBS})
add_executable(lcfit-test EXCLUDE_FROM_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/test_lcfit_cpp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lcfit_cpp.cc
  )
target_link_libraries(lcfit-test
  gsl
  gslcblas
  lcfit-static)