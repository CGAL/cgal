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

set(thingi10k_BLACKLIST_WITHOUT_MERGE_FACETS
    106789.stl
    # At point 4.0163683382116631 -2.094120689431076 0
    # There is a edge of the input mesh that is infinitely close to a vertex.
    )

function(CGAL_add_cdt3_test_from_Thingi10k data_name data_filename)
  set(options "ONLY_MERGE_FACETS")
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
  if(MY_ONLY_MERGE_FACETS)
    set(MY_ONLY_MERGE_FACETS ONLY_MERGE_FACETS)
  endif()
  CGAL_add_cdt3_from_off_test_aux(${data_name} ${THINGI10K_DATA_DIR} DATA_FILENAME ${data_filename}
     LABELS Thingi10K ${MY_LABELS}
     ${MY_TIMEOUT_KEYWORD} ${MY_TIMEOUT}
     ${MY_ONLY_MERGE_FACETS}
     )
endfunction()

foreach(thingi_file_name ${thingi10k_max_10k_solid})
  if(thingi_file_name IN_LIST thingi10k_BLACKLIST_WITHOUT_MERGE_FACETS)
    set(MY_ONLY_MERGE_FACETS ONLY_MERGE_FACETS)
  endif()
  get_filename_component(thingi_ID "${thingi_file_name}" NAME_WE)
  CGAL_add_cdt3_test_from_Thingi10k(Thingi10K_${thingi_ID} ${thingi_file_name}
      TIMEOUT 600 LABELS Thing10K_max_10k_solid ${MY_ONLY_MERGE_FACETS})
endforeach()
