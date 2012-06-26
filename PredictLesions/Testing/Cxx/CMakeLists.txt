
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

set(testname ${CLP}Test)
add_test(${testname} ${Launcher_Command} ${CMAKE_CURRENT_BINARY_DIR}/${CLP}Test
  --compare ${BASELINE}/outputLesionProb.nrrd ${TEMP}/outputLesionProb.nrrd
  ModuleEntryPoint
  --inputT1Volume ${INPUT}/t1_3mm.nrrd
  --inputT2Volume ${INPUT}/t2_3mm.nrrd
  --inputFLAIRVolume ${INPUT}/flair_3mm.nrrd
  --inputMaskVolume ${INPUT}/brain_mask_3mm.nrrd
  --inputT1RefVolume ${INPUT}/t1_3mm.nrrd
  --inputT2RefVolume ${INPUT}/t2_3mm.nrrd
  --inputFLAIRRefVolume ${INPUT}/flair_3mm.nrrd
  --inputMaskRefVolume ${INPUT}/brain_mask_3mm.nrrd
  --inputModel ${INPUT}/lesionSegmentation.model
  --inputLesionThreshold 10
  --outputLesionVolume ${TEMP}/outputLesionMask.nrrd
  --outputLesionProbVolume ${TEMP}/outputLesionProb.nrrd
  )
set_tests_properties ( ${testname} PROPERTIES TIMEOUT 100)
