## Add executable ##
add_executable(
  clamped_bar
  "main.cpp"
  ${SOURCES})

add_executable(
  optimize_clamped_bar
  "optimize.cpp"
  ${SRC_DIR}/clamped_bar_problem.cpp
  ${SOURCES})

message("Boost root dir: " ${BOOST_ROOT})


## Specify libreries that have to be linked to the bins ##
target_link_libraries(
  clamped_bar
  ${VTK_LIBRARIES}
  Qt5::Widgets
  # ${Boost_GRAPH_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  nlopt
  )

target_link_libraries(
  optimize_clamped_bar
  pagmo
  ${VTK_LIBRARIES}
  Qt5::Widgets
  ${Boost_GRAPH_LIBRARY}
  ${Boost_THREAD_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  nlopt
)

## Set binary output location  ##
# Do nothing atm

## Set target compiler options ##
set(DBG_FLAGS "-Wall")
set(REL_FLAGS "-O3" "-DNDEBUG")
set(COMPILER_FLAGS ${REL_FLAGS} "-std=c++14" "-pthread")
message("--CLAMPED BAR-- \nCompiler flags: " ${COMPILER_FLAGS} "\n")
target_compile_options(clamped_bar PUBLIC ${COMPILER_FLAGS} )
target_compile_options(optimize_clamped_bar PUBLIC ${COMPILER_FLAGS} )


