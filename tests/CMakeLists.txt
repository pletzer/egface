include_directories(${CMAKE_SOURCE_DIR})

set(LD_LIBRARY_PATH_NAME "LD_LIBRARY_PATH")
if (APPLE)
   set(LD_LIBRARY_PATH_NAME "DYLD_LIBRARY_PATH")
endif()

add_executable(testPointIntersectorLineMarginal
               "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLineMarginal.cxx")
target_link_libraries(testPointIntersectorLineMarginal
                      egface ${VTK_LIBRARIES})

add_executable(testSimplex
               "${CMAKE_SOURCE_DIR}/tests/testSimplex.cxx")
target_link_libraries(testSimplex
                      egface ${VTK_LIBRARIES})

add_executable(testSimplex2
               "${CMAKE_SOURCE_DIR}/tests/testSimplex2.cxx")
target_link_libraries(testSimplex2
                      egface ${VTK_LIBRARIES})

add_executable(testNodalBasis
               "${CMAKE_SOURCE_DIR}/tests/testNodalBasis.cxx")
target_link_libraries(testNodalBasis
                      egface ${VTK_LIBRARIES})

add_executable(testEdgeBasis
               "${CMAKE_SOURCE_DIR}/tests/testEdgeBasis.cxx")
target_link_libraries(testEdgeBasis
                      egface ${VTK_LIBRARIES})

add_executable(testFaceBasis
               "${CMAKE_SOURCE_DIR}/tests/testFaceBasis.cxx")
target_link_libraries(testFaceBasis
                      egface ${VTK_LIBRARIES})

add_executable(testSimplexIntersect
               "${CMAKE_SOURCE_DIR}/tests/testSimplexIntersect.cxx")
target_link_libraries(testSimplexIntersect
                      egface ${VTK_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(testSimplexIter
               "${CMAKE_SOURCE_DIR}/tests/testSimplexIter.cxx")
target_link_libraries(testSimplexIter
                      egface ${VTK_LIBRARIES})

message(STATUS "VTK_DIR=${VTK_DIR}")
#
# class tests
#
add_test(NAME testGrid
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testGrid.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugrid2d.vtk"
         --output "testGrid.vtk")
set_tests_properties(testGrid 
   PROPERTIES 
   ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../")

add_test(NAME testCellLocator
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testCellLocator.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/rbc_001.vtk")
set_tests_properties(testCellLocator
   PROPERTIES 
   ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../")

add_test(NAME testPointIntersectorLineMarginalCxx
         COMMAND "${CMAKE_BINARY_DIR}/tests/testPointIntersectorLineMarginal"
         "${CMAKE_SOURCE_DIR}/data/ugridex.vtk")
set_tests_properties(testPointIntersectorLineMarginalCxx
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "has [^013456789] points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )

add_test(NAME testPointIntersectorPoint0
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorPoint.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/cube.vtk"
         --p0 "0.5, 0.2, 0.3"
         )
set_tests_properties(testPointIntersectorPoint0
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "cell  0  has  1 intersection points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )

add_test(NAME testPointIntersectorPoint1
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorPoint.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/cube.vtk"
         --p0 "0.5, 0.8, 0.3"
         )
set_tests_properties(testPointIntersectorPoint1
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "cell  1  has  1 intersection points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )

add_test(NAME testPointIntersectorPoint2
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorPoint.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/cube.vtk"
         --p0 "0.5, 0.5, 0.7"
         )
set_tests_properties(testPointIntersectorPoint2 
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "cell  4  has  1 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  5  has  1 intersection points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../")

add_test(NAME testPointIntersectorLineX
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/cube.vtk"
         --p0 "0, 0.4, 0.2" --p1 "1, 0.4, 0.2"
         )
set_tests_properties(testPointIntersectorLineX 
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "cell  0  has  2 intersection points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )

add_test(NAME testPointIntersectorLineY
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/cube.vtk"
         --p0 "0.5, 0, 0.2" --p1 "0.5, 1, 0.2"
         )
set_tests_properties(testPointIntersectorLineY
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "cell  0  has  2 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  1  has  2 intersection points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )

add_test(NAME testPointIntersectorLineZ
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/cube.vtk"
         --p0 "0.5, 0.4, 0" --p1 "0.5, 0.4, 1"
         )
set_tests_properties(testPointIntersectorLineZ
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "cell  0  has  2 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  2  has  2 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  4  has  2 intersection points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )


add_test(NAME testPointIntersectorLineNoIntersection
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "-1, -1, 0" --p1 "-2, -3, 0"
         )
set_tests_properties(testPointIntersectorLineNoIntersection
   PROPERTIES 
   ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../")


add_test(NAME testPointIntersectorLine2
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "1, 0.5, -1" --p1 "1, 0.5, 7"
         )
set_tests_properties(testPointIntersectorLine2
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "has [^013456789] points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )

add_test(NAME testPointIntersectorLineMarginal
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "1, 0.5, 0" --p1 "1, 0.5, 6"
         )
set_tests_properties(testPointIntersectorLineMarginal
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "has [^013456789] points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )

add_test(NAME testPointIntersectorLineInside
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLine.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "1, 0.5, 0.001" --p1 "1, 0.5, 5.999"
         )
set_tests_properties(testPointIntersectorLineInside
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "has [^013456789] points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )

add_test(NAME testPointIntersectorTriangleSimple
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorTriangle.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/cube.vtk"
         --p0 "0.5, 0.2, 0.3" --p1 "0.5, 0.8, 0.3" --p2 "0.5, 0.5, 0.7"
         )
set_tests_properties(testPointIntersectorTriangleSimple
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "cell  0  has  4 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  1  has  4 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  2  has  4 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  3  has  4 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  4  has  3 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  5  has  3 intersection points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )

add_test(NAME testPointIntersectorTetrahedronSimple
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorTetrahedron.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/cube.vtk"
         --p0 "0, 0, 0" --p1 "1, 0.5, 0" --p2 "0, 1, 0" --p3 "0, 0.5, 1"
         )
set_tests_properties(testPointIntersectorTetrahedronSimple
                     PROPERTIES 
                     PASS_REGULAR_EXPRESSION "cell  0  has  6 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  1  has  6 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  2  has  9 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  3  has  9 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  4  has  8 intersection points"
                     PASS_REGULAR_EXPRESSION "cell  5  has  8 intersection points"
                     ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../"
                     )

# tests
add_test(NAME testPointIntersectorTriangleOutside
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorTriangle.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "-1, -1, 0" --p1 "-2, -1, 0" --p2 "-2, -3, 0"
         )
set_tests_properties(testPointIntersectorTriangleOutside
   PROPERTIES 
   ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../")

add_test(NAME testPointIntersectorTriangleInside
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorTriangle.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "1, 0.2, 0.1" --p1 "1, 0.2, 5" --p2 "1, 0.7, 3"
         )
set_tests_properties(testPointIntersectorTriangleInside
   PROPERTIES 
   ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../")

add_test(NAME testPointIntersectorTriangleMarginal
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorTriangle.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "1, 0, 0" --p1 "1, 0, 6" --p2 "1, 1, 3"
         )
set_tests_properties(testPointIntersectorTriangleMarginal
   PROPERTIES 
   ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../")

add_test(NAME testPointIntersectorTetrahedronOutside
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorTetrahedron.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "-1, -1, 0" --p1 "-2, -1, 0" --p2 "-2, -3, 0" --p3 "-4, -3, -1"
         )
set_tests_properties(testPointIntersectorTetrahedronOutside
   PROPERTIES 
   ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../")

add_test(NAME testPointIntersectorTetrahedronInside
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorTetrahedron.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "1, 0.2, 0.1" --p1 "1, 0.2, 5" --p2 "1, 0.7, 3" --p3 "1.5, 0.7, 3"
         )
set_tests_properties(testPointIntersectorTetrahedronInside
   PROPERTIES 
   ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../")

add_test(NAME testPointIntersectorTetrahedronMarginal
         COMMAND "${PYTHON_EXECUTABLE}" 
         "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorTetrahedron.py"
         --build_dir ${CMAKE_BINARY_DIR}
         --input "${CMAKE_SOURCE_DIR}/data/ugridex.vtk"
         --p0 "0, 0, 0" --p1 "0, 0, 6" --p2 "0, 1, 3" --p3 "2, 1, 3"
         )
set_tests_properties(testPointIntersectorTetrahedronMarginal
   PROPERTIES 
   ENVIRONMENT "${LD_LIBRARY_PATH_NAME}=${VTK_DIR}/../../")




