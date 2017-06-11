## Add executable ##
add_executable(thermo "main.cpp" ${SOURCES})

## Specify libraries that have to be linked to the binaries ##
target_link_libraries(
  thermo
  ${VTK_LIBRARIES}
  Qt5::Widgets
  ${Boost_GRAPH_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
)

## Set binary output location  ##
# Do nothing atm

## Set target compiler options ##
set(COMPILER_FLAGS "-std=c++14" "-Wall" "-g")
message("--TEST THERMO-- \nCompiler flags: " ${COMPILER_FLAGS} "\n")
target_compile_options(thermo PUBLIC ${COMPILER_FLAGS})
