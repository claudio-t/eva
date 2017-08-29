## Add executable ##
add_executable(frame2d "main.cpp" ${SOURCES})

## Specify libraries that have to be linked to the binaries ##
target_link_libraries(
  frame2d
  ${VTK_LIBRARIES}
  Qt5::Widgets
  ${Boost_GRAPH_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
)
  


## Set binary output location  ##
# Do nothing atm

## Set target compiler options ##
set(COMPILER_FLAGS "-std=c++14" "-Wall" "-g")
message("--TEST FRAME2D-- \nCompiler flags: " ${COMPILER_FLAGS} "\n")
target_compile_options(frame2d PUBLIC ${COMPILER_FLAGS})

