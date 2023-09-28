if(CGAL_CDT_3_NO_THINGI10K)
  return()
endif()
find_path(THINGI10K_DATA_DIR NAME 132423.stl
  HINTS ENV HOME
  PATH_SUFFIXES Downloads/Thingi10K/raw_meshes
  NO_DEFAULT_PATH
  NO_CMAKE_FIND_ROOT_PATH
  )
include(./thingi10k_max_10k_solid.cmake)

function(CGAL_add_cdt3_test_from_Thingi10k data_name data_filename)
  set(options "")
  set(oneValueArgs TIMEOUT)
  set(multiValueArgs LABELS)
  cmake_parse_arguments(PARSE_ARGV 2 "MY" "${options}" "${oneValueArgs}"
                      "${multiValueArgs}")
  if(MY_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unknown arguments specified: ${MY_UNPARSED_ARGUMENTS}")
  endif()
  if(MY_TIMEOUT)
    set(MY_TIMEOUT_KEYWORD TIMEOUT)
  endif()
  CGAL_add_cdt3_from_off_test_aux(${data_name} ${THINGI10K_DATA_DIR} DATA_FILENAME ${data_filename}
     LABELS Thingi10K ${MY_LABELS}
     ${MY_TIMEOUT_KEYWORD} ${MY_TIMEOUT}
     )
endfunction()

foreach(thingi_file_name ${thingi10k_max_10k_solid})
  get_filename_component(thingi_ID "${thingi_file_name}" NAME_WE)
  CGAL_add_cdt3_test_from_Thingi10k(Thingi10K_${thingi_ID} ${thingi_file_name} TIMEOUT 600 LABELS Thing10K_max_10k_solid)
endforeach()
