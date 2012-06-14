
foreach(dir flann
    flann/algorithms
    flann/io
    flann/mpi
    flann/nn
    flann/util
    flann/util/cuda
    )
  message("dir = ${dir}")
  file(GLOB include_files ${FLANN_SOURCE}/src/cpp/${dir}/*.h
    ${FLANN_SOURCE}/cpp/*.hpp)
  foreach(include_file ${include_files})
    get_filename_component(include_filename ${include_file} NAME)
    configure_file(${include_file}
      ${INCLUDE_TARGET}/${dir}/${include_filename} COPY_ONLY)
  endforeach()
endforeach()

