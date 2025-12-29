cmake_minimum_required (VERSION 3.13)
project(libMcPhase)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED 11)
cmake_policy(SET CMP0148 OLD)  # New policy gives error when cannot find shared libs

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

if (DEFINED ENV{CONDA_PREFIX})
    message(STATUS "Searching $ENV{CONDA_PREFIX} for Pybind11 and Eigen")
    file(GLOB pybind11_search $ENV{CONDA_PREFIX}/share/cmake/pybind11/pybind11Config.cmake)
    if (NOT pybind11_search STREQUAL "")
        cmake_path(GET pybind11_search PARENT_PATH pybind11_DIR)
    endif()
    file(GLOB eigen_search $ENV{CONDA_PREFIX}/share/eigen3/cmake/Eigen3Config.cmake)
    if (NOT eigen_search STREQUAL "")
        cmake_path(GET eigen_search PARENT_PATH Eigen3_DIR)
    endif()
endif()

if (EMSCRIPTEN)
    message(STATUS "Emscripten detected - setting Pybind to use cross-compiler")
    set(PYBIND11_USE_CROSSCOMPILING TRUE)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(pybind11 CONFIG)

add_subdirectory(src)

pybind11_add_module(${LIBMCPHASE_PYTHON_MODULE} MODULE)
target_include_directories(${LIBMCPHASE_PYTHON_MODULE} PUBLIC src/include)
add_subdirectory(src/libmcphase)
