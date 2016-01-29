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

add_library(egface SHARED
            egfGrid.cpp)

# enable ctest
enable_testing()

#
# class tests
#
add_test(NAME testGrid
         COMMAND "${PYTHON_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/tests/testGrid.py")

