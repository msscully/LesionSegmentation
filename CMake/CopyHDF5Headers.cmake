file(GLOB c_includes ${HDF5_SOURCE}/src/*.h)
file(GLOB cxx_includes ${HDF5_SOURCE}/c++/src/*.h)
file(GLOB config_includes ${HDF5_BUILD}/*.h)
# make the include directory, including subdirectories
file(MAKE_DIRECTORY ${INCLUDE_TARGET}/hdf5/cpp)

foreach(include_file ${c_includes} ${config_includes})
  get_filename_component(include_filename ${include_file} NAME)
  configure_file(${include_file} ${INCLUDE_TARGET}/hdf5/${include_filename} COPYONLY)
endforeach()

foreach(include_file ${cxx_includes})
  get_filename_component(include_filename ${include_file} NAME)
  configure_file(${include_file} ${INCLUDE_TARGET}/hdf5/cpp/${include_filename} COPYONLY)
endforeach()
