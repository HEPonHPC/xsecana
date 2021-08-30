cmake_minimum_required(VERSION 3.14)

project(XSecAna)
include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(ROOT REQUIRED)
find_package(Eigen3 REQUIRED)
if(NOT CAFANA_INCLUDE_DIR)
       set(CAFANA_INCLUDE_DIR "$ENV{CAFANACORE_DIR}include")
endif()

message(STATUS "ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}")
message(STATUS "CAFANA_INCLUDE_DIR ${CAFANA_INCLUDE_DIR}")



add_subdirectory(test)
