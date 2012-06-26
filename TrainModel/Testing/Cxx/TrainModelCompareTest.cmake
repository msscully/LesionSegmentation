# arguments checking
if( NOT TEST_PROGRAM )
  message( FATAL_ERROR "Require TEST_PROGRAM to be defined" )
endif( NOT TEST_PROGRAM )
if( NOT TEST_COMPARE_PROGRAM )
  message( FATAL_ERROR "Require TEST_COMPARE_PROGRAM to be defined" )
endif( NOT TEST_COMPARE_PROGRAM )
if( NOT TEST_BASELINE )
  message( FATAL_ERROR "Require TEST_BASELINE to be defined" )
endif( NOT TEST_BASELINE )
if( NOT TEST_INPUT_DIR )
  message( FATAL_ERROR "Require TEST_INPUT_DIR to be defined" )
endif( NOT TEST_INPUT_DIR )
if( NOT TEST_TEMP_OUTPUT )
  message( FATAL_ERROR "Require TEST_TEMP_OUTPUT to be defined" )
endif( NOT TEST_TEMP_OUTPUT )

# Run the compare program to make sure it built correctly
execute_process(
  COMMAND ${TEST_COMPARE_PROGRAM} --help
  RESULT_VARIABLE TEST_RESULT
  )

# if the return value is !=0 bail out
if( TEST_RESULT )
  message( FATAL_ERROR "Failed: Test compare program ${TEST_COMPARE_PROGRAM} won't run.\n${TEST_ERROR}" )
endif( TEST_RESULT )

# Check to see if the image we are comparing against exists.  We do this here to avoid a lengthy test for no reason.
if(NOT EXISTS ${TEST_BASELINE})
  message( FATAL_ERROR "Failed: Baseline file ${TEST_BASELINE} does not exist!\n")
endif( NOT EXISTS ${TEST_BASELINE})

# run the test program, capture the stdout/stderr and the result var
execute_process(
	COMMAND ${TEST_PROGRAM} ModuleEntryPoint --inputT1Volumes ${TEST_INPUT_DIR}/t1_3mm.nrrd --inputT2Volumes ${TEST_INPUT_DIR}/t2_3mm.nrrd --inputFLAIRVolumes ${TEST_INPUT_DIR}/flair_3mm.nrrd --inputMaskVolumes ${TEST_INPUT_DIR}/brain_mask_3mm.nrrd --inputLesionVolumes ${TEST_INPUT_DIR}/lesion_3mm.nrrd --inputPercentNonLesion 1 --outputModel ${TEST_TEMP_OUTPUT}/lesionSegmentation.model
  ERROR_VARIABLE TEST_ERROR
  RESULT_VARIABLE TEST_RESULT
  )

# if the return value is !=0 bail out
if( TEST_RESULT )
  message( FATAL_ERROR "Failed: Test program ${TEST_PROGRAM} exited != 0.\n${TEST_ERROR}" )
endif( TEST_RESULT )

# now compare the output with the reference
execute_process(
  COMMAND ${TEST_COMPARE_PROGRAM} --inputModel2 ${TEST_TEMP_OUTPUT}/lesionSegmentation.model --inputModel1 ${TEST_BASELINE}
  RESULT_VARIABLE TEST_RESULT
  )

# again, if return value is !=0 scream and shout
if( TEST_RESULT )
  message( FATAL_ERROR "Failed: The output of ${TEST_PROGRAM} did not match ${TEST_BASELINE}")
endif( TEST_RESULT )

# everything went fine...
message( "Passed: The output of ${TEST_PROGRAM} matches ${TEST_BASELINE}" )

