# egface

cmake_minimum_required (VERSION 2.8)
project(egface)

# need python
find_package(PythonInterp REQUIRED)
# need VTK include files and libraries
find_package(VTK REQUIRED)

# make sure the numpy python module can be imported
execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import numpy" ERROR_VARIABLE ERROR_NUMPY)
if (NOT ERROR_NUMPY STREQUAL "")
  message(FATAL_ERROR "You must have python-numpy installed!")
endif()

include(${VTK_USE_FILE})
add_library(egface SHARED
            egfGrid.cpp
            egfCellLocator.cpp
            egfPointIntersector.cpp)
target_link_libraries(egface ${VTK_LIBRARIES})

# enable ctest
enable_testing()

#
# class tests
#
add_test(NAME testGrid
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testGrid.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugrid2d.vtk"
         --output "testGrid.vtk")

add_test(NAME testCellLocator
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testCellLocator.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/rbc_001.vtk")

add_test(NAME testPointIntersectorLineNoIntersection
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "-1, -1, 0" --p1 "-2, -3, 0"
         )

add_test(NAME testPointIntersectorLine2
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "1, 0.5, -1" --p1 "1, 0.5, 7"
         )

add_test(NAME testPointIntersectorLineMarginal
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "1, 0.5, 0" --p1 "1, 0.5, 6"
         )

add_test(NAME testPointIntersectorLineInside
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "1, 0.5, 0.001" --p1 "1, 0.5, 5.999"
         )

