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

execute_process(
  COMMAND ${CMD}
  ERROR_VARIABLE err
  OUTPUT_VARIABLE output
  RESULT_VARIABLE error_result)

file(WRITE ${CERR} "${err}")

if(error_result)
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
