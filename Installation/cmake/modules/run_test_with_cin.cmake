if(NOT CMD OR (NOT ANDROID AND NOT CGAL_RUN_TESTS_THROUGH_SSH AND NOT EXISTS ${CMD}))
  message(FATAL_ERROR
    "The variable `CMD` should be defined to the test program to run!")
endif()

if(NOT CIN OR NOT EXISTS ${CIN})
  message(FATAL_ERROR
    "The variable `CIN` should be defined to the input file for the test!")
endif()

if(ANDROID)
  execute_process(
    COMMAND adb shell "cd ${ANDROID_DIR_PREFIX}${PROJECT_NAME} && ${ANDROID_DIR_PREFIX}${PROJECT_NAME}/${CMD}"
    INPUT_FILE ${CIN}
    RESULT_VARIABLE error_result)
elseif(CGAL_RUN_TESTS_THROUGH_SSH)
  execute_process(
    COMMAND ${ssh_executable} ${SSH_HOST} "cd ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME} && LD_LIBRARY_PATH=${CGAL_REMOTE_TEST_DIR_PREFIX}/lib ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}/${CMD}"
    INPUT_FILE ${CIN}
    RESULT_VARIABLE error_result)
else()
  execute_process(
    COMMAND ${CMD}
    INPUT_FILE ${CIN}
    RESULT_VARIABLE error_result)
endif()

if(error_result)
  message(SEND_ERROR
    "The test `${CMD} < ${CIN}` ended with the error code ${error_result}.")
endif()
