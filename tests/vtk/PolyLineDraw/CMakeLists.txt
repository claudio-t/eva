cmake_minimum_required(VERSION 2.8)

PROJECT(test_vtk)

find_package(VTK REQUIRED)
find_package(Qt5 COMPONENTS Widgets)
# find_package(Qt5Core REQUIRED)
# find_package(Qt5 REQUIRED)

include(${VTK_USE_FILE})

add_executable(test_vtk test_vtk)

foreach(l ${VTK_LIBRARIES})
  MESSAGE(${l})
endforeach()

# MESSAGE( STATUS "Qt5 LIBS:         " ${QT_LIBRARIES} )

# if(VTK_LIBRARIES)
target_link_libraries(test_vtk ${VTK_LIBRARIES} Qt5::Widgets)
# else()
  # target_link_libraries(test_vtk vtkHybrid vtkWidgets)
# endif()