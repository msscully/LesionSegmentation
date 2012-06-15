# build Flann as an external project
ExternalProject_Add(FLANN
  URL "http://people.cs.ubc.ca/~mariusm/uploads/FLANN/flann-1.7.1-src.zip"
  URL_MD5 d780795f523eabda7c7ea09c6f5cf235
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BINARY_DIR FLANN-build
  SOURCE_DIR FLANN
  CMAKE_ARGS
  ${CMAKE_COMMON_ARGS}
  -DBUILD_PYTHON_BINDINGS:BOOL=OFF
  -DBUILD_MATLAB_BINDINGS:BOOL=OFF
  -DBUILD_CUDA_LIB:BOOL=OFF
  DEPENDS HDF5
  )

set(INCLUDE_TARGET ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

#
# run script to copy the include files
ExternalProject_Add_Step(FLANN CopyHeaders
  COMMENT "Copying FLANN Headers"
  DEPENDEES build
  COMMAND ${CMAKE_COMMAND}
  -DINCLUDE_TARGET=${INCLUDE_TARGET}
  -DFLANN_BUILD=<BINARY_DIR>
  -DFLANN_SOURCE=<SOURCE_DIR>
  -P ${CMAKE_CURRENT_LIST_DIR}/CopyFLANNHeaders.cmake)

#
# use import_libraries (in ExtProjectSetup.cmake)
# to import the libraries and creat the library var
set(FLANN_LibNames flann_cpp_s flann_cpp_s-gd flann_s)
set(FLANN_LibDir ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

import_libraries(EXTPROJECT FLANN LIBNAMES ${FLANN_LibNames}
  LIBDIR ${FLANN_LibDir} LIBVARNAME FLANN_LIBRARIES)

#
# since the headers are sprinkled in several directories, loop
# to add them all
set(FLANN_IncludeDirs
  flann/algorithms flann/io  flan/mpi flann/nn  flan/tbb  flan/util
  flann/util/cuda flan)
include_directories(${INCLUDE_TARGET})
foreach(incdir ${FLANN_IncludeDirs})
  include_directories(${INCLUDE_TARGET}/${incdir})
endforeach()

set(FLANN_INCLUDE_DIR ${INCLUDE_TARGET})
IF(BUILD_SHARED_LIBS)
set(FLANN_LIBRARY $<TARGET_FILE:flann_cpp>)
ELSE(BUILD_SHARED_LIBS)
set(FLANN_LIBRARY $<TARGET_FILE:flann_cpp_s>)
ENDIF(BUILD_SHARED_LIBS)
