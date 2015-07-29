cmake_minimum_required (VERSION 2.6) 
project (clustering)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

find_package(Eigen2 REQUIRED)
include_directories(${Eigen2_INCLUDE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/../lib)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

link_directories(${CMAKE_SOURCE_DIR}/../build/lib/)

add_executable(example example.cpp)
add_definitions(-std=c++11)
target_link_libraries(example clustering)