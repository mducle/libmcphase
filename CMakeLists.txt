cmake_minimum_required (VERSION 3.13)
project(libMcPhase)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED 11)

set(LIBMCPHASE_PYTHON_MODULE libmcphase)

# The latest VS needs an extra library VCRUNTIME140_1.dll which is not in 
# the anaconda redistributable module (yet) - so disable it.
# https://developercommunity.visualstudio.com/content/problem/852548/vcruntime140-1dll-is-missing.html
if(WIN32)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -d2FH4-")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -d2FH4-")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -d2:-FH4-")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -d2:-FH4-")
endif(WIN32)

include(Eigen)
include_directories(${EIGEN3_INCLUDE_DIR})

include(PyBind11)
add_subdirectory(${PyBind11_DIR})

add_subdirectory(src)

pybind11_add_module(${LIBMCPHASE_PYTHON_MODULE} MODULE)
target_include_directories(${LIBMCPHASE_PYTHON_MODULE} PUBLIC src/include)
add_subdirectory(src/libmcphase)
