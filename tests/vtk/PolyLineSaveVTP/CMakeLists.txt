cmake_minimum_required(VERSION 2.8)

PROJECT(SaveToVTP)

find_package(VTK REQUIRED)
find_package(Qt5 COMPONENTS Widgets)

include(${VTK_USE_FILE})

add_executable(save_to_vtp save_to_vtp)

target_link_libraries(save_to_vtp ${VTK_LIBRARIES} Qt5::Widgets)
# target_link_libraries(save_to_vtp vtkHybrid vtkWidgets Qt5::Widgets)


