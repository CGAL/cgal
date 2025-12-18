if(CGAL_CDT_3_NO_THINGI10K)
  return()
endif()
find_path(THINGI10K_DATA_DIR NAME 132423.stl
  HINTS ENV HOME
  PATH_SUFFIXES Downloads/Thingi10K/raw_meshes
  NO_DEFAULT_PATH
  NO_CMAKE_FIND_ROOT_PATH
  REQUIRED
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

set(thingi10k_FAILED_WITH_SEGFAULT_CTest_20251002
1439534.stl
196123.stl
200695.stl
135777.stl
285604.stl
822697.stl
)

set(thingi10k_FAILED_CTest_20251002
100606.stl
100644.stl
101955.stl
109130.stl
116873.stl
116876.stl
135777.stl
139737.stl
1439534.stl
145329.stl
145330.stl
1505036.stl
1514900.stl
196121.stl
196122.stl
196123.stl
196126.stl
196127.stl
199814.stl
199818.stl
200695.stl
215991.stl
230152.stl
230153.stl
239188.stl
276937.stl
285604.stl
285605.stl
288352.stl
288353.stl
288354.stl
288355.stl
39182.stl
39245.stl
472050.stl
55278.stl
61418.stl
622000.stl
669962.stl
67817.stl
702204.stl
723893.stl
822697.stl
904476.stl
91474.stl
95796.stl
95797.stl
97515.stl
)

set(thingi10k_FAILED_WITH_MERGE_FACETS_CTest_20251028
139765.stl
1452677.stl
1452678.stl
1452679.stl
145329.stl
145330.stl
145331.stl
1505036.stl
1514900.stl
153100.stl
1652975.stl
1652976.stl
1706457.stl
186546.stl
196121.stl
196122.stl
196123.stl
196126.stl
196127.stl
196194.stl
199814.stl
199818.stl
206318.stl
215991.stl
230152.stl
230153.stl
237632.stl
239188.stl
255657.stl
255658.stl
276937.stl
285603.stl
286161.stl
288352.stl
288446.stl
360073.stl
362398.stl
37743.stl
383022.stl
39182.stl
39245.stl
39495.stl
39499.stl
40841.stl
41521.stl
42040.stl
44025.stl
44064.stl
44901.stl
472050.stl
50659.stl
51797.stl
57811.stl
61418.stl
61431.stl
622000.stl
62592.stl
62593.stl
65144.stl
65395.stl
65402.stl
669962.stl
68255.stl
702204.stl
70381.stl
71461.stl
723893.stl
72419.stl
726665.stl
77343.stl
84624.stl
90225.stl
906183.stl
91147.stl
91474.stl
93702.stl
93703.stl
95796.stl
95797.stl
97515.stl
97590.stl
97593.stl
99895.stl
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
  if(thingi_file_name IN_LIST thingi10k_FAILED_CTest_20251002)
    list(APPEND LABELS "CTest_20251002_failed")
  endif()
  if(thingi_file_name IN_LIST thingi10k_FAILED_WITH_SEGFAULT_CTest_20251002)
    list(APPEND LABELS "CTest_20251002_failed_segfault")
  endif()
  if(thingi_file_name IN_LIST thingi10k_FAILED_WITH_MERGE_FACETS_CTest_20251028)
    list(APPEND LABELS "CTest_20251028_failed_merge_facets")
  endif()
  get_filename_component(thingi_ID "${thingi_file_name}" NAME_WE)
  CGAL_add_cdt3_test_from_Thingi10k(Thingi10K_${thingi_ID} ${thingi_file_name}
      TIMEOUT 600 LABELS ${LABELS} ${MY_ONLY_MERGE_FACETS})
endforeach()
