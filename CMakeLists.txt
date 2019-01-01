cmake_minimum_required (VERSION 3.1)
project(libMcPhase)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_MACOSX_RPATH TRUE)

set(LIBMCPHASE_VERSION_MAJOR "0")
set(LIBMCPHASE_VERSION_MINOR "1")
set(LIBMCPHASE_VERSION_PATCH "0")
set(LIBMCPHASE_VERSION "${LIBMCPHASE_VERSION_MAJOR}.${LIBMCPHASE_VERSION_MINOR}.${LIBMCPHASE_VERSION_PATCH}")

set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED 11)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(Eigen)
include_directories(${EIGEN3_INCLUDE_DIR})

include(PyBind11)
add_subdirectory(${PyBind11_DIR} ${PyBind11_DIR})

include_directories(include)

add_library(singleion OBJECT singleion/cfpars.cpp singleion/cf1ion.cpp singleion/racah.cpp
                             singleion/ic_states.cpp singleion/ic_cfp.cpp singleion/ic_coulomb.cpp
                             singleion/ic_socf.cpp singleion/ic1ion.cpp)
set_property(TARGET singleion PROPERTY POSITION_INDEPENDENT_CODE ON)

add_subdirectory(pythonapi)

enable_testing()
add_test(NAME libMcPhasePython COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/bin 
         python3 -m unittest discover WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
