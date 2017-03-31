include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../lcfit_src
  ${CMAKE_CURRENT_SOURCE_DIR}/../lcfit_cpp_src)

add_executable(lcfit-test EXCLUDE_FROM_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/test_lcfit.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/test_lcfit_cpp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/test_lcfit2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/test_lcfit3.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cc)
target_link_libraries(lcfit-test
  lcfit_cpp-static)

add_test(NAME lcfit-test COMMAND lcfit-test)
