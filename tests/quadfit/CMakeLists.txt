add_executable (quadfit-sim quadfit.cc )
SET_TARGET_PROPERTIES(quadfit-sim
  PROPERTIES OUTPUT_NAME quadfit)
target_link_libraries (quadfit-sim cxxlsfit-lib)

add_test (quadfit ${CMAKE_CURRENT_SOURCE_DIR}/../../tools/basic_test.sh
  ${CMAKE_CURRENT_BINARY_DIR}/quadfit ${CMAKE_CURRENT_SOURCE_DIR})

