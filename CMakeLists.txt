cmake_minimum_required (VERSION 3.13)
project(libMcPhase)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED 11)

set(LIBMCPHASE_PYTHON_MODULE libmcphase)

include(Eigen)
include_directories(${EIGEN3_INCLUDE_DIR})

include(PyBind11)
add_subdirectory(${PyBind11_DIR})

add_subdirectory(src)

pybind11_add_module(${LIBMCPHASE_PYTHON_MODULE} MODULE src/libmcphase/libmcphase.cpp)

# Adds C++ libraries to Python module
target_include_directories(${LIBMCPHASE_PYTHON_MODULE} PUBLIC src/include)
target_link_libraries(${LIBMCPHASE_PYTHON_MODULE} PUBLIC ${LIBMCPHASE_CXX_LIBRARIES})
