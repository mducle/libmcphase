include(ExternalProject)

option(USE_SYSTEM_EIGEN "Use the system installed Eigen" OFF)

if(USE_SYSTEM_EIGEN)
  message(STATUS "Using system Eigen")
  find_package(Eigen3 REQUIRED)
else()
  message(STATUS "Using Eigen in ExternalProject")
  ExternalProject_Add(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/extern-eigen/download
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/extern-eigen/source
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/extern-eigen/install
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    TEST_COMMAND ""
  )
  set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/extern-eigen/source")
endif()
