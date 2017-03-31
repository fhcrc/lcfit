cmake_minimum_required (VERSION 2.8.1)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic")
set(CMAKE_C_FLAGS_DEBUG "-g -DLCFIT_DEBUG -DLCFIT4_VERBOSE -DLCFIT2_VERBOSE -DLCFIT3_VERBOSE -DLCFIT_AUTO_VERBOSE")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DVERBOSE")

add_subdirectory(lcfit_src)
add_subdirectory(lcfit_cpp_src)
add_subdirectory(test)
