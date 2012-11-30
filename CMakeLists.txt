cmake_minimum_required (VERSION 2.8.1)

project(lcfit CXX)

set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -pedantic")

# STS
add_subdirectory(sts)

set(LCFIT_COMPARE_CPP ${CMAKE_CURRENT_SOURCE_DIR}/src/compare.cc)
set(LCFIT_TEST_CPP ${CMAKE_CURRENT_SOURCE_DIR}/src/test.cc)
include_directories(src)
include_directories(sts/lib/include)

# LCFIT
add_executable(lcfit-compare ${LCFIT_COMPARE_CPP})
add_executable(lcfit-test ${LCFIT_TEST_CPP})
target_link_libraries(lcfit-test sts-static ${STS_PHYLO_LIBS})
target_link_libraries(lcfit-compare sts-static ${STS_PHYLO_LIBS})
