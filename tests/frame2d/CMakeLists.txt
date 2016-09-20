## Add executable ##
add_executable(frame2d "frame2d.cpp" ${SOURCES})

## Specify libraries that have to be linked to the binaries ##
target_link_libraries(frame2d ${VTK_LIBRARIES} Qt5::Widgets ${Boost_GRAPH_LIBRARY})

## Set binary output location  ##
# Do nothing atm

## Set target compiler options ##
target_compile_options(frame2d PUBLIC "-std=c++14" "-O2" "-g")
