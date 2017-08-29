## Add executable ##
add_executable(
  optimize-clamped-bar
  "optimize.cpp"
  ${SRC_DIR}/thermo_frame_problem.cpp
  ${SRC_DIR}/thermo_frame_problem_utils.cpp
  ${SOURCES})
add_executable(
  post-process
  "post_process.cpp"
  ${SRC_DIR}/thermo_frame_problem.cpp
  ${SRC_DIR}/thermo_frame_problem_utils.cpp
  ${SOURCES})

## Specify libreries that have to be linked to the bins ##
target_link_libraries(
  optimize-clamped-bar
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

target_link_libraries(
  post-process
  pagmo
  ${VTK_LIBRARIES}
  Qt5::Widgets
  # ${Boost_GRAPH_LIBRARY}
  ${Boost_SERIALIZATION_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  nlopt
  )

## Set binary output location  ##
# Do nothing atm

## Set target compiler options ##
set(DBG_FLAGS "-Wall" "-g")
set(REL_FLAGS "-O3" "-DNDEBUG")
set(COMPILER_FLAGS ${REL_FLAGS} "-std=c++14" "-pthread")
message("--CLAMPED BAR-- \nCompiler flags: " ${COMPILER_FLAGS} "\n")
target_compile_options(optimize-clamped-bar PUBLIC ${COMPILER_FLAGS} )
target_compile_options(post-process PUBLIC ${COMPILER_FLAGS} )



