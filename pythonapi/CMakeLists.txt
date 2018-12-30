message(STATUS "PYTHON_SITE_PACKAGES = ${PYTHON_SITE_PACKAGES_DIR}" )

pybind11_add_module(libMcPhase libmcphase.cpp $<TARGET_OBJECTS:singleion>)

set_target_properties (libMcPhase PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/libMcPhase)
add_custom_command (TARGET libMcPhase POST_BUILD
                    COMMAND ${CMAKE_COMMAND} ARGS -E make_directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/libMcPhase)
add_custom_command (TARGET libMcPhase POST_BUILD
                    COMMAND ${CMAKE_COMMAND} ARGS -E copy_if_different 
                    ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/libMcPhase/)
