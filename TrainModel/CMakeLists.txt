#
# CMake code common to both the SuperBuild
cmake_minimum_required(VERSION 2.8.7)

set(MODULE_NAME TrainModel)
set(MODULE_TITLE ${MODULE_NAME})
string(TOUPPER ${MODULE_NAME} MODULE_NAME_UPPER)

if(NOT Slicer_SOURCE_DIR)
  set(MODULE_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/4.1/Modules/TrainModel")
  set(MODULE_CATEGORY "Segmentation")
  set(MODULE_STATUS "Alpha")
  set(MODULE_CONTRIBUTORS "Mark Scully (BioMedical Mining, LLC)")
  set(MODULE_DESCRIPTION "This CLI is used to create new segmentation models based on a number of subjects'preprocessed T1, T2, FLAIR, and brain mask volumes as well as expert lesion segmentations.")
  set(MODULE_ACKNOWLEDGEMENTS "This work was supported by Partners and NAMIC.")
  set(MODULE_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/CLIExtensionTemplate/CLIExtensionTemplate.png?revision=19437&view=co")
  set(MODULE_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/4/4f/TrainModelInterface.png")
  set(MODULE_MAJOR_VERSION 0)
  set(MODULE_MINOR_VERSION 3)
  set(MODULE_PATCH_VERSION 0)
endif()

set(cli_module_include_directories
  ../nanoflann
  ../
  )

set(cli_module_SRCS
  ../LesionSegmentationModel.h
  ../nanoflann/nanoflann.hpp
  )

set(cli_module_target_libraries
  ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Configure XML description
configure_file(
  ${MODULE_NAME}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.xml
  @ONLY
  )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  #LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
  TARGET_LIBRARIES ${cli_module_target_libraries}
  INCLUDE_DIRECTORIES ${cli_module_include_directories}
  ADDITIONAL_SRCS ${cli_module_SRCS}
  #EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()
