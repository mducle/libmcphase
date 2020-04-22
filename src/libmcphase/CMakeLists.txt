
set(LIBMCPHASE_PYTHON_SOURCES
    libmcphase.cpp   # This must be first
    pycfpars.cpp
    pycf1ion.cpp
    pyic1ion.cpp
)
target_sources(${LIBMCPHASE_PYTHON_MODULE} PRIVATE ${LIBMCPHASE_PYTHON_SOURCES})
target_link_libraries(${LIBMCPHASE_PYTHON_MODULE} PUBLIC ${LIBMCPHASE_CXX_LIBRARIES})

