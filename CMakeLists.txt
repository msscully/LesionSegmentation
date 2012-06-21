#
# CMake code common to both the SuperBuild
cmake_minimum_required(VERSION 2.8.7)

set(MODULE_NAME LesionSegmentation)

if(NOT Slicer_SOURCE_DIR)
  set(MODULE_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Snapshot/Developers/Tutorials/ExtensionWriting")
  set(MODULE_CATEGORY "Segmentation")
  set(MODULE_STATUS "Alpha")
  set(MODULE_CONTRIBUTORS "Mark Scully (BioMedical Mining, LLC)")
  set(MODULE_DESCRIPTION "This is a set of CLI modules for creating lesion segmentation models and segmenting lesions using those models.")
  set(MODULE_ACKNOWLEDGEMENTS "This work was supported by Partners and NAMIC.")
  set(MODULE_MAJOR_VERSION 0)
  set(MODULE_MINOR_VERSION 1)
  set(MODULE_PATCH_VERSION 0)
endif()

if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

set(cli_module_include_directories
  )

set(cli_module_SRCS
  LesionSegmentationModel.h
  nanoflann/nanoflann.hpp
  )

set(cli_module_target_libraries
  ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
foreach(proj TrainModel PredictLesions)
  # Configure XML description
  configure_file(
    ${proj}.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/${proj}.xml
    @ONLY
    )

  slicerMacroBuildCLI(
    NAME ${proj}
    #LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
    TARGET_LIBRARIES ${cli_module_target_libraries}
    INCLUDE_DIRECTORIES ${cli_module_include_directories}
    ADDITIONAL_SRCS ${cli_module_SRCS}
    #EXECUTABLE_ONLY
    )
endforeach()

#-----------------------------------------------------------------------------
#if(BUILD_TESTING)
#  add_subdirectory(Testing)
#endif()

#-----------------------------------------------------------------------------

if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()
