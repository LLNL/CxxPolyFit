add_executable (quad2fit-sim quad2fit.cc )
SET_TARGET_PROPERTIES(quad2fit-sim
  PROPERTIES OUTPUT_NAME quad2fit)
target_link_libraries (quad2fit-sim cxxlsfit-lib)

add_test (quad2fit ${CMAKE_CURRENT_SOURCE_DIR}/../../tools/basic_test.sh
  ${CMAKE_CURRENT_BINARY_DIR}/quad2fit ${CMAKE_CURRENT_SOURCE_DIR})

