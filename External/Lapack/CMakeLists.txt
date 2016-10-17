file(GLOB lapack_SRC
  "*.f"
)
add_library(lapack ${LIBRARY_TYPE} ${lapack_SRC})

install (TARGETS lapack
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
