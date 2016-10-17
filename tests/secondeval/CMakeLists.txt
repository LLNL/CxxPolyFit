add_executable (secondeval-sim secondeval.cc )
SET_TARGET_PROPERTIES(secondeval-sim
  PROPERTIES OUTPUT_NAME secondeval)
target_link_libraries (secondeval-sim cxxlsfit-lib)

add_test (secondeval ${CMAKE_CURRENT_SOURCE_DIR}/../../tools/basic_test.sh
  ${CMAKE_CURRENT_BINARY_DIR}/secondeval ${CMAKE_CURRENT_SOURCE_DIR})

