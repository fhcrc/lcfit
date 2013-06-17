cmake_minimum_required (VERSION 2.8.1)

if(NOT NO_LCFIT)
add_definitions(-DUSE_LCFIT -DLCFIT_DEBUG)
add_subdirectory(lcfit_src)
add_subdirectory(lcfit_cpp_src)
endif()

project(phyml C)
set(CMAKE_C_FLAGS "-ansi -pedantic -Wall -std=c99")
set(CMAKE_C_FLAGS_RELEASE
  "${CMAKE_CFLAGS_RELEASE} -msse -fomit-frame-pointer -funroll-loops")
set(PHYML_SRC phyml_src/src)
set(PHYML_SOURCES
  ${PHYML_SRC}/utilities.c
  ${PHYML_SRC}/optimiz.c
  ${PHYML_SRC}/lk.c
  ${PHYML_SRC}/bionj.c
  ${PHYML_SRC}/models.c
  ${PHYML_SRC}/free.c
  ${PHYML_SRC}/help.c
  ${PHYML_SRC}/simu.c
  ${PHYML_SRC}/eigen.c
  ${PHYML_SRC}/pars.c
  ${PHYML_SRC}/alrt.c
  ${PHYML_SRC}/interface.c
  ${PHYML_SRC}/cl.c
  ${PHYML_SRC}/spr.c
  ${PHYML_SRC}/draw.c
  ${PHYML_SRC}/stats.c
  ${PHYML_SRC}/rates.c
  ${PHYML_SRC}/mcmc.c
  ${PHYML_SRC}/times.c
  ${PHYML_SRC}/tiporder.c
  ${PHYML_SRC}/mg.c
  ${PHYML_SRC}/m4.c
  ${PHYML_SRC}/io.c
  ${PHYML_SRC}/make.c
  ${PHYML_SRC}/nexus.c
  ${PHYML_SRC}/init.c
  ${PHYML_SRC}/xml.c
  ${PHYML_SRC}/sergeii.c
  ${PHYML_SRC}/mixt.c
  ${PHYML_SRC}/main.c
)
include_directories(${PHYML_SRC})

add_executable(phyml
  ${PHYML_SOURCES})

if(NOT NO_LCFIT)
target_link_libraries(phyml m lcfit-static gsl gslcblas)
include_directories(lcfit_src)
endif()
