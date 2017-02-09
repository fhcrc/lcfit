cmake_minimum_required (VERSION 2.8.1)

# Comparison
project(liblcfit_cpp C CXX)
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DVERBOSE")

set(LCFIT_LIB_CPP_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/gsl.h
  ${CMAKE_CURRENT_SOURCE_DIR}/lcfit_cpp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/lcfit_rejection_sampler.h)
set(LCFIT_LIB_CPP_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/gsl.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lcfit_cpp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/lcfit_rejection_sampler.cc)

add_library(lcfit_cpp-static STATIC ${LCFIT_LIB_CPP_FILES})
add_library(lcfit_cpp SHARED ${LCFIT_LIB_CPP_FILES})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lcfit_src)

target_link_libraries(lcfit_cpp-static
  ${GSL_LIBRARIES}
  lcfit-static
  nlopt)

target_link_libraries(lcfit_cpp
  ${GSL_LIBRARIES}
  lcfit
  nlopt)

# applications

set(LCFIT_COMPARE_CPP
  ${CMAKE_CURRENT_SOURCE_DIR}/lcfit_compare.cc)

set(LCFIT_APP_LIBS
  bpp-core
  bpp-seq
  bpp-phyl
  lcfit_cpp-static)

add_executable(lcfit-compare EXCLUDE_FROM_ALL ${LCFIT_COMPARE_CPP})
target_link_libraries(lcfit-compare ${LCFIT_APP_LIBS})
