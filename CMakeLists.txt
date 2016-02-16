# egface

cmake_minimum_required (VERSION 2.8)
project(egface)

# need python
find_package(PythonInterp REQUIRED)
# set VTK_DIR to the top directory
find_package(VTK REQUIRED)
# Lapack
find_package(LAPACK REQUIRED)

# make sure the numpy python module can be imported
execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import numpy" ERROR_VARIABLE ERROR_NUMPY)
if (NOT ERROR_NUMPY STREQUAL "")
  message(FATAL_ERROR "You must have python-numpy installed!")
endif()

include(${VTK_USE_FILE})
add_library(egface SHARED
            egfGrid.cpp
            egfCellLocator.cpp
            egfPointIntersector.cpp
            egfField.cpp
            Factorial.cpp
            Simplex.cpp
            SimplexIter.cpp
            SimplexIntersect.cpp
            NodalBasis.cpp
            EdgeBasis.cpp
            FaceBasis.cpp
            )

if(VTK_LIBRARIES)
  target_link_libraries(egface ${VTK_LIBRARIES} ${LAPACK_LIBRARIES})
else()
  target_link_libraries(egface vtkHybrid ${LAPACK_LIBRARIES})
endif()

# enable ctest
include(CTest)
enable_testing()

add_subdirectory(tests)
