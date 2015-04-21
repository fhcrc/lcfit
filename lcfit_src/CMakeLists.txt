cmake_minimum_required (VERSION 2.8)

# Build the LCFIT C library
project(libclfit C)
set(CMAKE_C_FLAGS "-std=c99 -Wall -pedantic")

set(LCFIT_LIB_C_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/lcfit.h
    ${CMAKE_CURRENT_SOURCE_DIR}/lcfit_select.h)
set(LCFIT_LIB_C_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/lcfit.c
    ${CMAKE_CURRENT_SOURCE_DIR}/lcfit_select.c)
add_library(lcfit-static STATIC ${LCFIT_LIB_C_FILES})
add_library(lcfit SHARED ${LCFIT_LIB_C_FILES})
target_link_libraries(lcfit-static
  gsl
  gslcblas
  nlopt)
target_link_libraries(lcfit
  gsl
  gslcblas
  nlopt)

install(TARGETS lcfit DESTINATION lib)
install(FILES ${LCFIT_LIB_C_HEADERS} DESTINATION include)
