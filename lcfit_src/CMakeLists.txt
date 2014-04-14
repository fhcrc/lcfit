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

set(LCFIT_LINK_LIBRARIES m gsl gslcblas)
target_link_libraries(lcfit-static ${LCFIT_LINK_LIBRARIES} nlopt)
target_link_libraries(lcfit ${LCFIT_LINK_LIBRARIES} nlopt-shared)

install(TARGETS lcfit DESTINATION lib)
install(FILES ${LCFIT_LIB_C_HEADERS} DESTINATION include)
