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
    #106789.stl
    # At point 4.0163683382116631 -2.094120689431076 0
    # There is a edge of the input mesh that is infinitely close to a vertex.
    )

set(thingi10k_FAILED_WITH_MERGE_FACETS
    1053875.stl
    112926.stl
    1147177.stl
    118295.stl
    123787.stl
    126284.stl
    135777.stl
    1452672.stl
    1514900.stl
    1514901.stl
    1514903.stl
    1514904.stl
    162100.stl
    162336.stl
    1743322.stl
    186554.stl
    186559.stl
    237737.stl
    239182.stl
    255657.stl
    285604.stl
    285605.stl
    288353.stl
    288354.stl
    288355.stl
    383022.stl
    43663.stl
    45550.stl
    464846.stl
    472002.stl
    472004.stl
    472050.stl
    472190.stl
    500116.stl
    520645.stl
    54467.stl
    55278.stl
    57811.stl
    57812.stl
    59229.stl
    67817.stl
    69260.stl
    71920.stl
    74780.stl
    90224.stl
    904476.stl
    904480.stl
    96046.stl
    96457.stl
    96659.stl
    97503.stl
)

set(thingi10k_FAILED_WITH_MERGE_FACETS_CTest_20240222_2201
40119.stl
40985.stl
41360.stl
44903.stl
47732.stl
55262.stl
55278.stl
57812.stl
58439.stl
67817.stl
77342.stl
80084.stl
90224.stl
92118.stl
97503.stl
112926.stl
135777.stl
162336.stl
186544.stl
186559.stl
225958.stl
285604.stl
285605.stl
288353.stl
288354.stl
288355.stl
375273.stl
442387.stl
464846.stl
904476.stl
904480.stl
1053875.stl
1452672.stl
1505023.stl
1514904.stl
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
     LABELS ${MY_LABELS}
     ${MY_TIMEOUT_KEYWORD} ${MY_TIMEOUT}
     ${MY_ONLY_MERGE_FACETS}
     )
endfunction()

foreach(thingi_file_name ${thingi10k_max_10k_solid})

  if(thingi_file_name IN_LIST thingi10k_BLACKLIST_WITHOUT_MERGE_FACETS)
    set(MY_ONLY_MERGE_FACETS ONLY_MERGE_FACETS)

    unset(MY_ONLY_MERGE_FACETS)
  endif()
  set(LABELS Thingi10K Thingi10K_max_10k_solid)
  if(thingi_file_name IN_LIST thingi10k_FAILED_WITH_MERGE_FACETS)
    list(APPEND LABELS "Thingi10K_FAIL")
  endif()
  if(thingi_file_name IN_LIST thingi10k_FAILED_WITH_MERGE_FACETS_CTest_20240222_2201)
    list(APPEND LABELS "CTest_20240222_2201_failed_merge_facets")
  endif()
  get_filename_component(thingi_ID "${thingi_file_name}" NAME_WE)
  CGAL_add_cdt3_test_from_Thingi10k(Thingi10K_${thingi_ID} ${thingi_file_name}
      TIMEOUT 600 LABELS ${LABELS} ${MY_ONLY_MERGE_FACETS})
endforeach()
