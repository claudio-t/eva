cmake_minimum_required(VERSION 2.8)
 
project(EVA)
set(EVA_MAIN_DIR "/home/claudio/eva")


## Find required packages ##
find_package(VTK  REQUIRED)
find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Boost COMPONENTS graph REQUIRED)

## Run cmake within vtk dir ##
include(${VTK_USE_FILE})

## Include Dirs ##
include_directories("/home/claudio/eva/include")

## Specify target bin and source files ##
add_executable(eva src/main.cpp)

## Specify libreries that have to be linked to the bin ##
target_link_libraries(eva ${VTK_LIBRARIES} Qt5::Widgets ${Boost_GRAPH_LIBRARY})

## Select clang as c++ compiler ##
set(CMAKE_CXX_COMPILER "/usr/bin/clang++") 

## Set compiler options ##
target_compile_options(eva PUBLIC "-std=c++14" "-Wall")