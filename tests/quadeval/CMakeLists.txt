add_executable (quadeval-sim quadeval.cc )
SET_TARGET_PROPERTIES(quadeval-sim
  PROPERTIES OUTPUT_NAME quadeval)
target_link_libraries (quadeval-sim cxxlsfit-lib)

add_test (quadeval ${CMAKE_CURRENT_SOURCE_DIR}/../../tools/basic_test.sh
  ${CMAKE_CURRENT_BINARY_DIR}/quadeval ${CMAKE_CURRENT_SOURCE_DIR})

