if(NOT DEFINED CMD)
#  message("CMAKE_ARGC: ${CMAKE_ARGC}")
#  message("CMAKE_ARGV0: ${CMAKE_ARGV0}")
#  message("CMAKE_ARGV1: ${CMAKE_ARGV1}")
#  message("CMAKE_ARGV2: ${CMAKE_ARGV2}")
#  message("CMAKE_ARGV3: ${CMAKE_ARGV3}")
#  message("CMAKE_ARGV4: ${CMAKE_ARGV4}")
  foreach(n RANGE 4 ${CMAKE_ARGC})
    list(APPEND CMD ${CMAKE_ARGV${n}})
  endforeach()
endif()
#message("run_cmd_redirection, the CMD list is: ${CMD}")
if(NOT CERR)
  message(FATAL_ERROR
    "The variable `CERR` should be defined to the output error file!")
endif()

# Create the file before using it
file(WRITE ${CERR})

# Execute the command ${CMD} with stderr redirected to the file ${CERR}
execute_process(
  COMMAND ${CMD}
  ERROR_FILE "${CERR}"
  OUTPUT_VARIABLE output
  RESULT_VARIABLE error_result)

if(error_result)
  if(CMD2)
    file(REMOVE ${CERR})
    execute_process(COMMAND ${CMD2})
    message(SEND_ERROR)
  else()
    file(READ ${CERR} err)
    file(REMOVE ${CERR})
    string(REPLACE ";" " " CMD_STR "${CMD}")
    message(SEND_ERROR
"The command
  ${CMD_STR} > ${CERR}
ended with the error code ${error_result},
the following output:
${output}
and the following error output:
${err}"
)
  endif()
endif()
