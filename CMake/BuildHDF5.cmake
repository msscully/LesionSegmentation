# build hdf5 as an external project
ExternalProject_Add(HDF5
  URL "http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8.7/src/hdf5-1.8.7.tar.gz"
  URL_MD5 37711d4bcb72997e93d495f97c76c33a  
  INSTALL_COMMAND ""
  UPDATE_COMMAND ""
  BINARY_DIR HDF5-build
  SOURCE_DIR HDF5
  CMAKE_ARGS
  ${CMAKE_COMMON_ARGS}
  -DHDF5_BUILD_CPP_LIB:BOOL=TRUE
  -DHDF5_BUILD_HL_LIB:BOOL=TRUE
#  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
  )

# after the build, copy the files to CMAKE_CURRENT_BINARY_DIR/include
set(INCLUDE_TARGET ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

ExternalProject_Add_Step(HDF5 CopyHeaders
  COMMENT "Copying HDF5 Headers"
  DEPENDEES build
  COMMAND ${CMAKE_COMMAND}
  -DINCLUDE_TARGET=${INCLUDE_TARGET}
  -DHDF5_BUILD=<BINARY_DIR>
  -DHDF5_SOURCE=<SOURCE_DIR>
  -P ${CMAKE_CURRENT_LIST_DIR}/CopyHDF5Headers.cmake)


#
# use import_libraries (in ExtProjectSetup.cmake)
# to import the libraries and creat the library var
set(HDF5_LibNames hdf5 hdf5_cpp hdf5_hl_cpp)
set(HDF5_LibDir ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

import_libraries(EXTPROJECT HDF5 LIBNAMES ${HDF5_LibNames}
  LIBDIR ${HDF5_LibDir} LIBVARNAME HDF5_LIBRARIES)

#
# The include directories
include_directories(${INCLUDE_TARGET}/hdf5
  ${INCLUDE_TARGET}/hdf5/cpp)

