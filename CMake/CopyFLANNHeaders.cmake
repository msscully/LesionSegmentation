
foreach(dir flann
    flan/tbb
    flann/algorithms
    flann/io
    flann/mpi
    flann/nn
    flann/util
    flann/util/cuda
    )
  message("dir = ${dir}")
  file(GLOB include_files ${FLANN_SOURCE}/src/cpp/${dir}/*.h
                          ${FLANN_SOURCE}/src/cpp/${dir}/*.hpp)
  foreach(include_file ${include_files})
    get_filename_component(include_filename ${include_file} NAME)
    set(destination ${INCLUDE_TARGET}/${dir}/${include_filename})

    message("copying ${include_file} to ${destination}")
    configure_file(${include_file} ${destination} COPY_ONLY)
  endforeach()
endforeach()

