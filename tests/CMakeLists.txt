include_directories(${CMAKE_SOURCE_DIR})
add_executable(testPointIntersectorLineMarginal
               "${CMAKE_SOURCE_DIR}/tests/testPointIntersectorLineMarginal.cxx")
target_link_libraries(testPointIntersectorLineMarginal
                      egface ${VTK_LIBRARIES})

