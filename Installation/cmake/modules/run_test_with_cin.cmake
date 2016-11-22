if(NOT CMD OR NOT EXISTS ${CMD})
  message(FATAL_ERROR
    "The variable `CMD` should be defined to the test program to run!")
endif()

if(NOT CIN OR NOT EXISTS ${CIN})
  message(FATAL_ERROR
    "The variable `CIN` should be defined to the input file for the test!")
endif()

execute_process(
  COMMAND ${CMD}
  INPUT_FILE ${CIN}
  RESULT_VARIABLE error_result)

if(error_result)
  message(SEND_ERROR
    "The test `${CMD} < ${CIN}` ended with the error code ${error_result}.")
endif()
