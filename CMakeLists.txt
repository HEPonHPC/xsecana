cmake_minimum_required(VERSION 3.14)

project(XSecAna)
include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(ROOT "6.18.00" REQUIRED)
if(NOT EIGEN3_INCLUDE_DIR)
       find_package(Eigen3 "3.4.0" REQUIRED)
endif()
set(EIGEN3_UNSUPPORTED_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR}/../unsupported/)

if(NOT CAFANA_INCLUDE_DIR)
       set(CAFANA_INCLUDE_DIR "$ENV{CAFANACORE_DIR}include")
endif()

message(STATUS "ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}")
message(STATUS "CAFANA_INCLUDE_DIR ${CAFANA_INCLUDE_DIR}")
message(STATUS "EIGEN3_INCLUDE_DIR ${EIGEN3_INCLUDE_DIR}")

set(CMAKE_CXX_FLAGS "-Wno-unknown-warning-option")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

include_directories(
        ./include/
        ${ROOT_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
)

add_subdirectory(test)
add_subdirectory(src)
