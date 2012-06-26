
#-----------------------------------------------------------------------------
set(Launcher_Command ${Slicer_LAUNCH_COMMAND})
set(BASELINE ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/Baseline)
set(INPUT ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/Input)
set(TEMP "${CMAKE_BINARY_DIR}/Testing/Temporary")

set(CLP ${MODULE_NAME})

#-----------------------------------------------------------------------------
add_executable(${CLP}Test ${CLP}Test.cxx)
target_link_libraries(${CLP}Test ${CLP}Lib)
set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

set(CLP CompareLesionModel)
set(${CLP}_SOURCE ${CLP}.cxx ${CMAKE_SOURCE_DIR}/LesionSegmentationModel.h)
generateclp(${CLP}_SOURCE ${CLP})
add_executable(${CLP} ${${CLP}_SOURCE})
set(CompareEXE ${CMAKE_CURRENT_BINARY_DIR}/${CLP})


add_test(NAME TrainModelTest COMMAND ${CMAKE_COMMAND}
	-D TEST_PROGRAM=${CMAKE_CURRENT_BINARY_DIR}/TrainModelTest
        -D TEST_COMPARE_PROGRAM=${CompareEXE}
	-D TEST_BASELINE=${INPUT}/lesionSegmentation.model
	-D TEST_INPUT_DIR=${INPUT}
        -D TEST_TEMP_OUTPUT=${TEMP}
	-P ${CMAKE_CURRENT_SOURCE_DIR}/TrainModelCompareTest.cmake
  )
set_tests_properties ( ${testname} PROPERTIES TIMEOUT 100)
