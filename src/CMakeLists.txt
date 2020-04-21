
include_directories(include)

set(SINGLEION_SOURCES
    singleion/cfpars.cpp
    singleion/cf1ion.cpp
    singleion/racah.cpp
    singleion/ic_states.cpp
    singleion/ic_cfp.cpp
    singleion/ic_coulomb.cpp
    singleion/ic_socf.cpp
    singleion/ic_tensorops.cpp
    singleion/ic1ion.cpp
)

# Adds the source modules as object libraries
add_library(singleion OBJECT ${SINGLEION_SOURCES})
set_property(TARGET singleion PROPERTY POSITION_INDEPENDENT_CODE ON)

# Put all the source modules into a single interface library for linking with Pybind11
set(LIBMCPHASE_CXX_LIBRARIES 
    singleion
CACHE INTERNAL "")

#add_subdirectory(pythonapi)
#enable_testing()
#add_test(NAME libMcPhasePython COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}/bin 
#         python3 -m unittest discover WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
