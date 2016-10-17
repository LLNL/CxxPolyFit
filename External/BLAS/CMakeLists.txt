file(GLOB blas_SRC
  "*.f"
)
add_library(blas ${LIBRARY_TYPE} ${blas_SRC})

install (TARGETS blas
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
