## Add executable ##
add_executable(eva "main.cpp" ${SOURCES})

## Specify libraries that have to be linked to the binaries ##
target_link_libraries(eva ${VTK_LIBRARIES} Qt5::Widgets ${Boost_GRAPH_LIBRARY})

## Set binary output location  ##
# Do nothing atm

## Set target compiler options ##
target_compile_options(eva PUBLIC "-std=c++14" "-Wall")
