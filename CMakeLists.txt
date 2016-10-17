cmake_minimum_required (VERSION 2.8)
project (cxxlsfit)

# enable testing
enable_testing ()

set (PROJECT_NAME cxxlsfit)
# The version number.
set (cxxlsfit_VERSION_MAJOR 0)
set (cxxlsfit_VERSION_MINOR 1)
set (cxxlsfit_VERSION_PATCH 0)

set (PROJECT_INCLUDE_DIR "${PROJECT_BINARY_DIR}/include")
set (PROJECT_BIN_DIR "${PROJECT_BINARY_DIR}/bin}")
set (PROJECT_LIB_DIR "${PROJECT_BINARY_DIR}/lib}")

set (EXTRA_LD_FLAGS "")


# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/Src/cxxlsfitCmakeConfig.h.in"
  "${PROJECT_INCLUDE_DIR}/cxxlsfitCmakeConfig.h"
  )

#########################################################
# Setup SYSTEM type, useful for some #defines in CXXLSFIT
#########################################################

if("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
add_definitions(-DLINUX)

elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
add_definitions(-DWINDOWS)

elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "MacOSX")
add_definitions(-DMACOS)

elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "cygwin")
add_definitions(-DCYGWIN)

elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "IRIX")
add_definitions(-DIRIX)

elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "HP-UX")
add_definitions(-DHPUX)

elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "AIX")
add_definitions(-DAIX)

elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "OFS/1")
add_definitions(-DALPHA)

elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "SunOS")
add_definitions(-DSUNOS)

else("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")

add_definitions(-DUSER)
endif("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")



#########################################################
# Setup the RPATH Junk
#########################################################
# use, i.e. don't skip the full RPATH for the build tree
if(BUILD_SHARED)

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_SKIP_RPATH FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

endif(BUILD_SHARED)
#########################################################
#END RPATH SETUP
#########################################################

#########################################################
# Special build types and extra flags
#########################################################
#SET(CMAKE_BUILD_TYPE nosse)
SET(CMAKE_CXX_FLAGS_NOSSE "-O0 -g -mno-sse")
SET(CMAKE_C_FLAGS_NOSSE "-O0 -g -mno-sse")
SET(CMAKE_Fortran_FLAGS_NOSSE "-O0 -g -mno-sse")

#SET(CMAKE_BUILD_TYPE testing)
SET(CMAKE_CXX_FLAGS_TESTING "-mfpmath=sse -msse2")
SET(CMAKE_C_FLAGS_TESTING "-mfpmath=sse -msse2")
SET(CMAKE_Fortran_FLAGS_TESTING "-mfpmath=sse -msse2")

#SET(CMAKE_BUILD_TYPE debug)
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -mfpmath=sse -msse2")
SET(CMAKE_C_FLAGS_DEBUG "-O0 -g -mfpmath=sse -msse2")
SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -mfpmath=sse -msse2")
 
SET(CMAKE_CXX_FLAGS_RELEASE "-O2")
SET(CMAKE_C_FLAGS_RELEASE "-O2")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")


# add the binary tree to the search path for include files
# so that we will find cxxlsfitCmakeConfig.h
include_directories("${PROJECT_INCLUDE_DIR}")
include_directories("Src")

# add the executable
file(GLOB cxxlsfit_SRC
  "Src/*.cc"
)

file(GLOB cxxlsfit_HDRS
  "Src/*.hh"
)

###############################################
#SETUP THE MAIN CXXLSFIT LIBRARY AND EXECUTABLES
###############################################
option (BUILD_SHARED "Build all the libraries as shared libs, and build libcxxlsfit.so" ON)
if (BUILD_SHARED)
  set(LIBRARY_TYPE SHARED)
else (BUILD_SHARED)
  set (LIBRARY_TYPE STATIC)
  SET(CMAKE_EXE_LINKER_FLAGS -static)
endif(BUILD_SHARED)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY )
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")


add_library (cxxlsfit-lib ${LIBRARY_TYPE} ${cxxlsfit_SRC} ${cxxlsfit_HDRS})
add_executable (cxxlsfit-bin Src/Main/cxxlsfit.cc)

set (CXXLSFIT_NAME "cxxlsfit")

if(BUILD_SHARED)

SET_TARGET_PROPERTIES(cxxlsfit-lib
  PROPERTIES
  OUTPUT_NAME ${CXXLSFIT_NAME})

SET_TARGET_PROPERTIES(cxxlsfit-bin
  PROPERTIES 
  OUTPUT_NAME ${CXXLSFIT_NAME}
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

else(BUILD_SHARED)

SET_TARGET_PROPERTIES(cxxlsfit-lib
  PROPERTIES
  OUTPUT_NAME ${CXXLSFIT_NAME}
  LINK_FLAGS "-static-libgcc  -static-libgfortran")

SET_TARGET_PROPERTIES(cxxlsfit-bin
  PROPERTIES
  OUTPUT_NAME ${CXXLSFIT_NAME}
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  LINK_FLAGS "-static-libgcc  -static-libgfortran")

endif(BUILD_SHARED)


###############################################
#SETUP 3RD PARTY LIBRARIES
###############################################


###################################################
# Lapack and BLAS are weird special cases, because we may find
# them or build it ourselves, so it gets this special section
###################################################
#Lapack and BLAS are required, but we should only build them
#if the aren't already installed  (Or if we're doing a release. Can't expect users to have them.)
#(Ordering is a little silly to get the link order right)


find_package(LAPACK)
if(LAPACK_FOUND AND NOT CMAKE_BUILD_TYPE MATCHES Release AND NOT MINGW)
  set (POLYFIT_EXTRA_LIBS ${POLYFIT_EXTRA_LIBS} ${LAPACK_LIBRARIES}) 
else()
  enable_language(Fortran)
  include_directories ("${PROJECT_SOURCE_DIR}/External/Lapack")
  add_subdirectory (External/Lapack) 
  set (POLYFIT_EXTRA_LIBS ${POLYFIT_EXTRA_LIBS} lapack)
endif()

find_package(BLAS)
if(BLAS_FOUND AND NOT CMAKE_BUILD_TYPE MATCHES Release AND NOT MINGW) #Build blas for windows and release
  set (POLYFIT_EXTRA_LIBS ${POLYFIT_EXTRA_LIBS} "${BLAS_LIBRARIES}" "gfortran") 
else()
  enable_language(Fortran)
  include_directories ("${PROJECT_SOURCE_DIR}/External/BLAS")
  add_subdirectory (External/BLAS) 
  set (POLYFIT_EXTRA_LIBS ${POLYFIT_EXTRA_LIBS} blas)
endif()

###################################################
# Link cxxlsfit library and the binary to the library
###################################################

target_link_libraries (cxxlsfit-lib ${POLYFIT_EXTRA_LIBS})
if("${CMAKE_SYSTEM_NAME}" STREQUAL "MacOSX")
  target_link_libraries (cxxlsfit-bin cxxlsfit-lib)  #I think OSX doesn't have "whole-archive?"  Linux and Windows need it
else("${CMAKE_SYSTEM_NAME}" STREQUAL "MacOSX")
  target_link_libraries (cxxlsfit-bin -Wl,-whole-archive cxxlsfit-lib  -Wl,-no-whole-archive)
endif("${CMAKE_SYSTEM_NAME}" STREQUAL "MacOSX")


###################################################
# Main cxxlsfit files installation section
###################################################
# add the install targets
install (TARGETS cxxlsfit-bin cxxlsfit-lib 
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)
install (FILES 
    "${PROJECT_BINARY_DIR}/include/cxxlsfitCmakeConfig.h"  
    ${cxxlsfit_HDRS}        
    DESTINATION include)

#######################################################
# With mingGw you have to include libgcc etc. in the package.  So do that here.
#######################################################
if( MINGW )
    get_filename_component( Mingw_Path ${CMAKE_CXX_COMPILER} PATH )
    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${Mingw_Path}/libgcc_s_dw2-1.dll ${Mingw_Path}/libstdc++-6.dll ${Mingw_Path}/libgfortran-3.dll ${Mingw_Path}/libquadmath-0.dll)
endif( MINGW )

include( InstallRequiredSystemLibraries )

if( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS )
    install( PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS} DESTINATION bin COMPONENT System )
endif( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS )

###################################################
# Packaging
###################################################
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_VERSION_MAJOR "${cxxlsfit_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${cxxlsfit_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${cxxlsfit_VERSION_PATCH}")
set (CPACK_PACKAGE_CONTACT       "tong10@llnl.gov")

if("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
  SET(CPACK_GENERATOR NSIS ZIP) #DEB and RPM are also options on appropriate systems
  SET(CPACK_NSIS_MODIFY_PATH ON)
  SET(CPACK_NSIS_CREATE_ICONS "CreateShortCut '\$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\${CXXLSFIT_NAME}.lnk' '\$INSTDIR\\\\bin\\\\${CXXLSFIT_NAME}.exe'") 

else()
  SET(CPACK_GENERATOR TGZ STGZ) #DEB and RPM are also options on appropriate systems
endif()

include (CPack)

add_subdirectory (tests/quadeval) 
add_subdirectory (tests/quadfit) 
add_subdirectory (tests/quad2eval) 
add_subdirectory (tests/quad2fit)  
add_subdirectory (tests/secondeval) 
add_subdirectory (tests/test3D) 

# enable dashboard scripting
include (CTest)
