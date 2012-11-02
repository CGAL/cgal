set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( NOT CGAL_VERSION_UTILS_FILE_INCLUDED )
  set( CGAL_VERSION_UTILS_FILE_INCLUDED 1 )

#  
# Given a version string of the form "major.[minor.[patch.[tweak]]]"
# decomposes it into components
#
macro( VERSION_DECOMPOSE v major minor patch tweak )

  string(REPLACE "." ";" VERSION_DECOMPOSE_LIST ${v} )

  list( LENGTH VERSION_DECOMPOSE_LIST VERSION_DECOMPOSE_LIST_LEN )
  
  if ( VERSION_DECOMPOSE_LIST_LEN GREATER 0 )
    list( GET VERSION_DECOMPOSE_LIST 0 ${major} )
  else()
    set ( ${major} -1 )
  endif()  
  
  if ( VERSION_DECOMPOSE_LIST_LEN GREATER 1 )
    list( GET VERSION_DECOMPOSE_LIST 1 ${minor} )
  else()
    set ( ${minor} -1 )
  endif()  
  
  if ( VERSION_DECOMPOSE_LIST_LEN GREATER 2 )
    list( GET VERSION_DECOMPOSE_LIST 2 ${patch} )
  else()
    set ( ${patch} -1 )
  endif()  
  
  if ( VERSION_DECOMPOSE_LIST_LEN GREATER 3 )
    list( GET VERSION_DECOMPOSE_LIST 3 ${tweak} )
  else()
    set ( ${tweak} -1 )
  endif()  
  
endmacro()

#  
# Given two version string of the form "major.[minor.[patch.[tweak]]]"
# returns TRUE if they are equal, FALSE otherwise.
#
macro( IS_VERSION_EQUAL a b r )

  VERSION_DECOMPOSE( ${a} _IVE_a_major _IVE_a_minor _IVE_a_patch _IVE_a_tweak )
  VERSION_DECOMPOSE( ${b} _IVE_b_major _IVE_b_minor _IVE_b_patch _IVE_b_tweak )
  
  set ( ${r} FALSE )  
  
  if ( _IVE_a_major EQUAL ${_IVE_b_major} )
    if ( _IVE_a_minor EQUAL ${_IVE_b_minor} )
      if ( _IVE_a_patch EQUAL ${_IVE_b_patch} )
        if ( _IVE_a_tweak EQUAL ${_IVE_b_tweak} )
          set ( ${r} TRUE )  
        endif()  
      endif()
    endif()
  endif()

endmacro()

#  
# Given two version string of the form "major.[minor.[patch.[tweak]]]"
# returns TRUE if the first is smaller than the second, FALSE otherwise.
#
macro( IS_VERSION_LESS a b r )

  VERSION_DECOMPOSE( ${a} _IVL_a_major _IVL_a_minor _IVL_a_patch _IVL_a_tweak )
  VERSION_DECOMPOSE( ${b} _IVL_b_major _IVL_b_minor _IVL_b_patch _IVL_b_tweak )
  
  set ( ${r} FALSE )  
  
  if ( _IVL_a_major LESS ${_IVL_b_major} )
    set ( ${r} TRUE )  
  elseif( _IVL_a_major EQUAL ${_IVL_b_major})  
    if ( _IVL_a_minor LESS ${_IVL_b_minor} )
      set ( ${r} TRUE )  
    elseif( _IVL_a_minor EQUAL ${_IVL_b_minor} )  
      if ( _IVL_a_patch LESS ${_IVL_b_patch} )
        set ( ${r} TRUE )  
      elseif( _IVL_a_patch EQUAL ${_IVL_b_patch} )   
        if ( _IVL_a_tweak LESS ${_IVL_b_tweak} )
          set ( ${r} TRUE )  
        endif()  
      endif()
    endif()
  endif()

endmacro()

#  
# Given two version string of the form "major.[minor.[patch.[tweak]]]"
# returns TRUE if the first is greater than the second, FALSE otherwise.
#
macro( IS_VERSION_GREATER a b r )

  IS_VERSION_LESS( ${a} ${b} _IVG_less )
  
  if ( _IVG_less )
    set( ${r} FALSE )
  else()
  
    IS_VERSION_EQUAL( ${a} ${b} _IVG_eq )
    
    if ( _IVG_eq )
      set( ${r} FALSE )
    else()
      set( ${r} TRUE )
    endif()
    
  endif()  
  
endmacro()














#
#                                    -= TESTING =-
#


macro ( TEST_VERSION_DECOMPOSE v expected_major expected_minor expected_patch expected_tweak )

  VERSION_DECOMPOSE( ${v} major minor patch tweak )
  
  set ( OK 0 )
  
  if ( major EQUAL "${expected_major}" )
    if ( minor EQUAL "${expected_minor}" )  
      if ( patch EQUAL "${expected_patch}" ) 
        if ( tweak EQUAL "${expected_tweak}" ) 
          set ( OK 1 )
        endif()  
      endif()  
    endif()  
  endif()  

  if ( OK )
    message( STATUS "correct - ${v} -> ${major}, ${minor}, ${patch}, ${tweak} " ) 
  else()
    message( STATUS "FAILED  - ${v} -> ${major}, ${minor}, ${patch}, ${tweak} " ) 
  endif()
  
endmacro()

macro ( TEST_VERSION_COMPARISON op v0 v1 expected )

  if ( "${op}" STREQUAL "<" )
    IS_VERSION_LESS( ${v0} ${v1} result )
  elseif ( "${op}" STREQUAL ">" )
    IS_VERSION_GREATER( ${v0} ${v1} result )
  else()
    IS_VERSION_EQUAL( ${v0} ${v1} result )
  endif()  
  
  if ( result STREQUAL ${expected} )
    message( STATUS "correct - ${v0} ${op} ${v1} => ${result}" )   
  else()
    message( STATUS "FAILED  - ${v0} ${op} ${v1} => ${result}" ) 
  endif()

endmacro()

if ( UNIT_TEST_VERSION_UTILS )

  TEST_VERSION_DECOMPOSE("1.2.3.4" 1 2 3 4 )
  TEST_VERSION_DECOMPOSE("1.2.3" 1 2 3 -1 )
  TEST_VERSION_DECOMPOSE("1.2" 1 2 -1 -1 )
  TEST_VERSION_DECOMPOSE("1" 1 -1 -1 -1 )
  
  TEST_VERSION_COMPARISON( "==" "1.2.3.4" "1.2.3.4" TRUE )
  TEST_VERSION_COMPARISON( "==" "1.2.3"   "1.2.3"   TRUE )
  TEST_VERSION_COMPARISON( "==" "1.2"     "1.2"     TRUE )
  TEST_VERSION_COMPARISON( "==" "1"       "1"       TRUE )
  
  TEST_VERSION_COMPARISON( "==" "1.2.3.4" "1.2.3"  FALSE )
  TEST_VERSION_COMPARISON( "==" "1.2.3"   "1.2"    FALSE )
  TEST_VERSION_COMPARISON( "==" "1.2"     "1"      FALSE )
  
  TEST_VERSION_COMPARISON( "==" "1.2.3.4" "1.2.3.5" FALSE )
  TEST_VERSION_COMPARISON( "==" "1.2.3"   "1.2.4"   FALSE )
  TEST_VERSION_COMPARISON( "==" "1.2"     "1.3"     FALSE )
  TEST_VERSION_COMPARISON( "==" "1"       "2"       FALSE )
  
  TEST_VERSION_COMPARISON( "<" "1.2.3.4" "1.2.3.4" FALSE )
  TEST_VERSION_COMPARISON( "<" "1.2.3"   "1.2.3"   FALSE )
  TEST_VERSION_COMPARISON( "<" "1.2"     "1.2"     FALSE )
  TEST_VERSION_COMPARISON( "<" "1"       "1"       FALSE )
  
  TEST_VERSION_COMPARISON( "<" "1.2.3.4" "1.2.3.5" TRUE )
  TEST_VERSION_COMPARISON( "<" "1.2.3.4" "1.2.4.5" TRUE )
  TEST_VERSION_COMPARISON( "<" "1.2.3.4" "1.3.4.5" TRUE )
  TEST_VERSION_COMPARISON( "<" "1.2.3.4" "2.3.4.5" TRUE )
  
  TEST_VERSION_COMPARISON( "<" "1.2.3.4" "1.2.4" TRUE )
  TEST_VERSION_COMPARISON( "<" "1.2.3.4" "1.3"   TRUE )
  TEST_VERSION_COMPARISON( "<" "1.2.3.4" "2"     TRUE )
  
  TEST_VERSION_COMPARISON( "<" "1.2.3" "1.2.4" TRUE )
  TEST_VERSION_COMPARISON( "<" "1.2"   "1.3"   TRUE )
  TEST_VERSION_COMPARISON( "<" "1"     "2"     TRUE )
  
  TEST_VERSION_COMPARISON( "<" "1.2.3.6" "1.2.4.5" TRUE )
  TEST_VERSION_COMPARISON( "<" "1.2.5.6" "1.3.4.5" TRUE )
  TEST_VERSION_COMPARISON( "<" "1.4.5.6" "2.3.4.5" TRUE )

  TEST_VERSION_COMPARISON( ">" "1.2.3.4" "1.2.3.4" FALSE )
  TEST_VERSION_COMPARISON( ">" "1.2.3"   "1.2.3"   FALSE )
  TEST_VERSION_COMPARISON( ">" "1.2"     "1.2"     FALSE )
  TEST_VERSION_COMPARISON( ">" "1"       "1"       FALSE )
  
  TEST_VERSION_COMPARISON( ">" "1.2.3.5" "1.2.3.4" TRUE )
  
endif()

endif()
