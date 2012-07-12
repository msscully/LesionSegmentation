cmake_minimum_required(VERSION 2.8.7)

set(EXTENSION_NAME LesionSegmentation)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/4.1/Extensions/LesionSegmentation")
  set(EXTENSION_CATEGORY "Segmentation")
  set(EXTENSION_STATUS "Alpha")
  set(EXTENSION_CONTRIBUTORS "Mark Scully (BioMedical Mining, LLC)")
  set(EXTENSION_DESCRIPTION "This is a set of CLI modules for creating lesion segmentation models and segmenting lesions using those models.")
  set(EXTENSION_ACKNOWLEDGEMENTS "This work was supported by Partners and NAMIC.")
  set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/CLIExtensionTemplate/CLIExtensionTemplate.png?revision=19437&view=co")
  set(EXTENSION_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/4/4f/TrainModelInterface.png http://www.slicer.org/slicerWiki/images/2/22/PredictModelInterface.png")
  set(EXTENSION_MAJOR_VERSION 0)
  set(EXTENSION_MINOR_VERSION 3)
  set(EXTENSION_PATCH_VERSION 0)
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
add_subdirectory(TrainModel)
add_subdirectory(PredictLesions)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()
