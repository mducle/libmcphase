include(ExternalProject)

option(USE_SYSTEM_EIGEN "Use the system installed Eigen" OFF)

if(USE_SYSTEM_EIGEN)
  message(STATUS "Using system Eigen")
  find_package(Eigen3 REQUIRED)
else()
  message(STATUS "Using Eigen in ExternalProject")

  # Download and unpack Eigen at configure time
  configure_file(${CMAKE_SOURCE_DIR}/cmake/Eigen.in ${CMAKE_BINARY_DIR}/extern-eigen/CMakeLists.txt)

  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/extern-eigen )
  execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/extern-eigen )

  set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/extern-eigen/source")
endif()
