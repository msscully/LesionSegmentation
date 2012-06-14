#
# CMake code common to both the SuperBuild
cmake_minimum_required(VERSION 2.8.7)

include(ExternalProject)
# get local files to include
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMake)

#
# the SuperBuild part of the project
if(NOT LesionSegmentation_NestedProject)
project(SuperBuild)
# Sets up the common arguments passed into ExternalProjects.
include(ExtProjectSetup)


#
# The problem with ExternalProjects is that if
# you need anything from the project BEFORE you
# configure your own project -- a program to run,
# a generated config file -- you have to build your
# own project as an External Build.
#
# In this case, you need GenerateCLP to convert
# the XML describing the command line parameters.
# So, we build ITK and SlicerExecutionModel
# then build the project.
include(BuildITK)
include(BuildSEM)

ExternalProject_add(LesionSegmentation
  DEPENDS ITKv4 SlicerExecutionModel
  DOWNLOAD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR LesionSegmentation-build
  CMAKE_ARGS
  ${CMAKE_COMMON_ARGS}
  -DLesionSegmentation_NestedProject:BOOL=ON
  -DITK_DIR:PATH=${ITK_DIR}
  -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
)

else(NOT LesionSegmentation_NestedProject)
# the actual project part of the build
project(LesionSegmentation)

# Sets up the common arguments passed into ExternalProjects.
include(ExtProjectSetup)


# build HDF5
include(BuildHDF5)
# build FLANN
include(BuildFLANN)

set(MODULE_NAME TrainModel)

if(NOT Slicer_SOURCE_DIR)
  set(MODULE_HOMEPAGE
    "http://www.slicer.org/slicerWiki/index.php/Documentation/Snapshot/Developers/Tutorials/ExtensionWriting")
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

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${SlicerExecutionModel_USE_FILE})
include(${SlicerExecutionModel_CMAKE_DIR}/SEMMacroBuildCLI.cmake)

set(SIMPLIFIED "ON")
if(SIMPLIFIED)

# Configure XML description
configure_file(
  ${MODULE_NAME}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.xml
  @ONLY
  )
set(cli_module_target_libraries
  ${ITK_LIBRARIES}
  ${FLANN_LIBRARIES}
  ${HDF5_LIBRARIES}
  )

SEMMacroBuildCLI(NAME PredictLesions
  ADDITIONAL_SRCS LesionSegmentationModel.h
  TARGET_LIBRARIES ${cli_modules_target_libraries}
  CLI_LIBRARY_WRAPPER_CXX CLIWrapper.cxx
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  INSTALL_RUNTIME_DESTINATION ${CMAKE_INSTALL_RUNTIME_DESTINATION}
  INSTALL_LIBRARY_DESTINATION ${CMAKE_CLI_INSTALL_LIBRARY_DESTINATION}
  INSTALL_ARCHIVE_DESTINATION ${CMAKE_CLI_INSTALL_ARCHIVE_DESTINATION}
)

else(SIMPLIFIED)

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
  ${FLANN_LIBRARIES}
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

set(MODULE_NAME PredictLesions)

if(NOT Slicer_SOURCE_DIR)
  set(MODULE_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Snapshot/Developers/Tutorials/ExtensionWriting")
  set(MODULE_CATEGORY "Segmentation")
  set(MODULE_STATUS "Alpha")
  set(MODULE_CONTRIBUTORS "Mark Scully (BioMedical Mining, LLC)")
  set(MODULE_DESCRIPTION "This is a CLI module for sementing lesions using a preexisting model.")
  #set(MODULE_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/CLIExtensionTemplate/CLIExtensionTemplate.png?revision=19437&view=co")
  #set(MODULE_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-r19441-CLIExtensionTemplate-screenshot.png http://wiki.slicer.org/slicerWiki/images/1/1e/Slicer-r19441-CLIExtensionTemplate-screenshot-2.png")
  set(MODULE_ACKNOWLEDGEMENTS "This work was supported by Partners and NAMIC.")
  set(MODULE_MAJOR_VERSION 0)
  set(MODULE_MINOR_VERSION 1)
  set(MODULE_PATCH_VERSION 0)
endif()

# Configure XML description
configure_file(
  ${MODULE_NAME}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.xml
  @ONLY
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

if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()
endif(SIMPLIFIED)
endif(NOT LesionSegmentation_NestedProject)
