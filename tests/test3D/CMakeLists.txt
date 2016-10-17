add_executable (test3D-sim test3D.cc )
SET_TARGET_PROPERTIES(test3D-sim
  PROPERTIES OUTPUT_NAME test3D)
target_link_libraries (test3D-sim cxxlsfit-lib)

add_test (test3D ${CMAKE_CURRENT_SOURCE_DIR}/../../tools/basic_test.sh
  ${CMAKE_CURRENT_BINARY_DIR}/test3D ${CMAKE_CURRENT_SOURCE_DIR})

