# used below in macro...
include(CMakeParseArguments)

macro(setifempty)
  if("${${ARGV0}}" STREQUAL "")
    set(${ARGV})
  endif()
endmacro(setifempty)
#-----------------------------------------------------------------------------
# the idea is to as much as possible to avoid
# installing projects to get them in the right place.
# setting this variables steers the ext projects to
# dump their libs and executables at the top level bin and lib
setifempty(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
setifempty(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
setifempty(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
setifempty(CMAKE_BUNDLE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
setifempty(CMAKE_INCLUDE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
# If you build pretty much anything on Linux without -fPIC, you lose
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(CMAKE_C_FLAGS "-fPIC ${CMAKE_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "-fPIC ${CMAKE_C_FLAGS}")
endif()

#
# establish a list of variables that always need passing to
# ext projects
set(CMAKE_COMMON_ARGS_LIST
  MAKECOMMAND:STRING
  CMAKE_SKIP_RPATH:BOOL
  CMAKE_BUILD_TYPE:STRING
  BUILD_SHARED_LIBS:BOOL
  CMAKE_CXX_COMPILER:PATH
  CMAKE_CXX_FLAGS_RELEASE:STRING
  CMAKE_CXX_FLAGS_DEBUG:STRING
  CMAKE_CXX_FLAGS:STRING
  CMAKE_C_COMPILER:PATH
  CMAKE_C_FLAGS_RELEASE:STRING
  CMAKE_C_FLAGS_DEBUG:STRING
  CMAKE_C_FLAGS:STRING
  CMAKE_SHARED_LINKER_FLAGS:STRING
  CMAKE_EXE_LINKER_FLAGS:STRING
  CMAKE_MODULE_LINKER_FLAGS:STRING
  CMAKE_GENERATOR:STRING
  CMAKE_EXTRA_GENERATOR:STRING
  CMAKE_INSTALL_PREFIX:PATH
  CMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH
  CMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH
  CMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH
  CMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH
  CMAKE_INCLUDE_OUTPUT_DIRECTORY:PATH
  CTEST_NEW_FORMAT:BOOL
  MEMORYCHECK_COMMAND_OPTIONS:STRING
  MEMORYCHECK_COMMAND:PATH
  CMAKE_SHARED_LINKER_FLAGS:STRING
  CMAKE_EXE_LINKER_FLAGS:STRING
  CMAKE_MODULE_LINKER_FLAGS:STRING
  SITE:STRING
  BUILDNAME:STRING
  )

#
# pick the list apart and build the arguments for
# ext project configuration
set(CMAKE_COMMON_ARGS)
foreach(arg ${CMAKE_COMMON_ARGS_LIST})
  string(REPLACE ":" ";" varname_and_vartype ${arg})
  set(target_info_list ${target_info_list})
  list(GET varname_and_vartype 0 _varname)
  list(GET varname_and_vartype 1 _vartype)
  list(APPEND CMAKE_COMMON_ARGS -D${_varname}:${_vartype}=${${_varname}})
endforeach()

# this is a first guess, I presume that
# on windows the library name is dependent on
# whether it's debug or release.
if(NOT WIN32 OR MINGW)
  set(lib_prefix lib)
  if(BUILD_SHARED_LIBS)
    set(lib_suffix .a)
  else()
    if(APPLE)
      set(lib_suffix .dylib)
    else()
      set(lib_suffix .so)
    endif()
  endif()
else(NOT WIN32)
  set(lib_prefix)
  if(BUILD_SHARED_LIBS)
    set(lib_suffix .dll)
  else()
    set(lib_suffix .lib)
  endif()
endif()

#
# handle importing libraries
# CMake doesn't know about libraries built outside
# the current project context. 'importing' libraries
# gives them first class targets you can use in
# this projects.
macro(import_libraries)
  set(oneValueArgs EXTPROJECT LIBDIR LIBVARNAME)
  set(multiValueArgs LIBNAMES)
  # parse the arguments
  cmake_parse_arguments(IMPLIBS "" # option args
    "${oneValueArgs}" # one value args
    "${multiValueArgs}" # multivalue args
    ${ARGN})

  # clear the Proj_LIBRARIES VAR
  set(${IMPLIBS_LIBVARNAME})
  #
  # import each lib
  foreach(lib ${IMPLIBS_LIBNAMES})
    # import as proper type of lib
    if(BUILD_SHARED_LIBS)
      add_library(${lib} SHARED IMPORTED)
    else()
      add_library(${lib} STATIC IMPORTED)
    endif()
    # add to Proj_LIBRARIES list
    list(APPEND ${IMPLIBS_LIBVARNAME} ${lib})
    # depend on the external project that
    # created this library
    add_dependencies(${lib}
      ${IMPLIBS_EXTPROJECT})
    # connect imported lib to actual file
    # created by EXTPROJECT
    set(libfilename
      ${IMPLIBS_LIBDIR}/${lib_prefix}${lib}${lib_suffix})
    set_property(TARGET ${lib} PROPERTY
      IMPORTED_LOCATION ${libfilename})
  endforeach()
endmacro()
