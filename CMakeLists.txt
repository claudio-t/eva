## Project config ##
cmake_minimum_required(VERSION 2.8)
project(EVA)
# set(CMAKE_VERBOSE_MAKEFILE ON)

## Setup aux vars ##
set(EVA_MAIN_DIR ${PROJECT_SOURCE_DIR})
set(SRC_DIR ${EVA_MAIN_DIR}/src)

## Select clang as c++ compiler ##
set(CMAKE_CXX_COMPILER "/usr/bin/clang++") 

## Force boost find module to use in usr/local/include-lib ##
set(BOOST_ROOT "/usr/local")
set(Boost_NO_SYSTEM_PATHS ON)

## Find required packages ##
find_package(VTK  REQUIRED)
find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Boost COMPONENTS graph thread serialization system program_options REQUIRED)

## Setup libraries ##
# set(EVA_LINK_LIBS
#   ${VTK_LIBRARIES}
#   Qt5::Widgets
#   ${Boost_GRAPH_LIBRARY})

## Setup VTK internal variables ##
include(${VTK_USE_FILE})

## Include Dirs ##
include_directories("include")

## Specify target bin and source files ##
set(SOURCES ${SRC_DIR}/vtk.cpp ${SRC_DIR}/utils.cpp ${SRC_DIR}/pagmo_utils.cpp)

## Specify subfolders ##
add_subdirectory("tests")
add_subdirectory("pagmo")


