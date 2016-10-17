add_executable (quad2eval-sim quad2eval.cc )
SET_TARGET_PROPERTIES(quad2eval-sim
  PROPERTIES OUTPUT_NAME quad2eval)
target_link_libraries (quad2eval-sim cxxlsfit-lib)

add_test (quad2eval ${CMAKE_CURRENT_SOURCE_DIR}/../../tools/basic_test.sh
  ${CMAKE_CURRENT_BINARY_DIR}/quad2eval ${CMAKE_CURRENT_SOURCE_DIR})

