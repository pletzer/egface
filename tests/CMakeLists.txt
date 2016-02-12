include_directories(${CMAKE_SOURCE_DIR})

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




