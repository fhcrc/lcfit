cmake_minimum_required (VERSION 2.8.1)

include_directories("nlopt-2.4.1/api")
add_subdirectory(nlopt-2.4.1)
add_subdirectory(lcfit_src)
add_subdirectory(lcfit_cpp_src)
