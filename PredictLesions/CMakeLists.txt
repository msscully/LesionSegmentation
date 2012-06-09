cmake_minimum_required(VERSION 2.8.7)

set(MODULE_NAME PredictLesions)

if(NOT Slicer_SOURCE_DIR)
  set(MODULE_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Snapshot/Developers/Tutorials/ExtensionWriting")
  set(MODULE_CATEGORY "Segmentation")
  set(MODULE_STATUS "Alpha")
  set(MODULE_CONTRIBUTORS "Mark Scully (BioMedical Mining, LLC)")
  set(MODULE_DESCRIPTION "This is a CLI module for creating lesion segmentation models.")
  #set(MODULE_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/CLIExtensionTemplate/CLIExtensionTemplate.png?revision=19437&view=co")
  #set(MODULE_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-r19441-CLIExtensionTemplate-screenshot.png http://wiki.slicer.org/slicerWiki/images/1/1e/Slicer-r19441-CLIExtensionTemplate-screenshot-2.png")
  set(MODULE_ACKNOWLEDGEMENTS "This work was supported by Partners and NAMIC.")
  set(MODULE_MAJOR_VERSION 0)
  set(MODULE_MINOR_VERSION 1)
  set(MODULE_PATCH_VERSION 0)
endif()

if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

# Configure XML description
configure_file(
  ${MODULE_NAME}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.xml
  @ONLY
  )

#-----------------------------------------------------------------------------
# FLANN (required)
find_package(Flann REQUIRED)

#-----------------------------------------------------------------------------

set(cli_module_include_directories
	${FLANN_INCLUDE_DIRS}
  )

set(cli_module_SRCS
	LesionSegmentationModel.h
  )

set(cli_module_target_libraries
  ${ITK_LIBRARIES}
  ${FLANN_LIBRARY}
  )

#-----------------------------------------------------------------------------
slicerMacroBuildCLI(
  NAME ${MODULE_NAME}
  #LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
  TARGET_LIBRARIES ${cli_module_target_libraries}
  INCLUDE_DIRECTORIES ${cli_module_include_directories}
  ADDITIONAL_SRCS ${cli_module_SRCS}
  EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
#if(BUILD_TESTING)
#  add_subdirectory(Testing)
#endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()
