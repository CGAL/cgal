#.rst:
# CGAL_SetupBoost
# ---------------
#
# The module searchs for the `Boost` headers and library, by calling
#
# .. code-block:: cmake
#
#    find_package(Boost)
#
# and defines the function :command:`use_CGAL_Boost_support`.

if ( CGAL_Boost_Setup )
  return()
endif()
set ( CGAL_Boost_Setup TRUE )

include(${CMAKE_CURRENT_LIST_DIR}/CGAL_TweakFindBoost.cmake)

function(CGAL_detect_if_Boost_Thread_is_required)
  get_property(PROPERTY_CGAL_requires_Boost_Thread_IS_SET
    GLOBAL PROPERTY CGAL_requires_Boost_Thread SET)
  if(PROPERTY_CGAL_requires_Boost_Thread_IS_SET)
    get_property(CGAL_requires_Boost_libs
      GLOBAL PROPERTY CGAL_requires_Boost_Thread)
  else()
    set ( CGAL_requires_Boost_libs TRUE )
    if ( DEFINED  MSVC_VERSION AND "${MSVC_VERSION}" GREATER 1800)
      set ( CGAL_requires_Boost_libs FALSE )
    else()
      try_run( CGAL_test_cpp_version_RUN_RES CGAL_test_cpp_version_COMPILE_RES
        "${CMAKE_BINARY_DIR}"
        "${CGAL_MODULES_DIR}/config/support/CGAL_test_cpp_version.cpp"
        RUN_OUTPUT_VARIABLE CGAL_cplusplus)
      message(STATUS "__cplusplus is ${CGAL_cplusplus}")
      if(NOT CGAL_test_cpp_version_RUN_RES)
        set ( CGAL_requires_Boost_libs FALSE )
        message(STATUS "  --> Do not link with Boost.Thread")
      endif()
    endif()
  endif()
  set_property(GLOBAL PROPERTY CGAL_requires_Boost_Thread ${CGAL_requires_Boost_libs})
  set(CGAL_requires_Boost_libs ${CGAL_requires_Boost_libs} PARENT_SCOPE)
endfunction(CGAL_detect_if_Boost_Thread_is_required)

CGAL_detect_if_Boost_Thread_is_required()
if (CGAL_requires_Boost_libs)
  find_package( Boost 1.48 REQUIRED thread system )
else()
  find_package( Boost 1.48 REQUIRED )
endif()

if(Boost_FOUND)
  if(DEFINED Boost_DIR AND NOT Boost_DIR)
    # Unset that cache variable that is set in the cache by FindBoost
    # (while it was searching for boost-cmake).
    unset(Boost_DIR CACHE)
    set(Boost_NO_BOOST_CMAKE TRUE CACHE INTERNAL "Avoid future search of boost-cmake")
  endif()
endif()

message( STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}" )
message( STATUS "Boost libraries:    ${Boost_LIBRARIES}" )

set ( CGAL_USE_BOOST 1 )


#.rst:
# Provided Functions
# ^^^^^^^^^^^^^^^^^^
#
# .. command:: use_CGAL_Boost_support
#
#    Link the target with the `Boost` libraries::
#
#      use_CGAL_Boost_support( target [INTERFACE] )
#
#    If the option ``INTERFACE`` is passed, the dependencies are
#    added using :command:`target_link_libraries` with the ``INTERFACE``
#    keyword, or ``PUBLIC`` otherwise.

function(use_CGAL_Boost_support target)
  if(ARGV1 STREQUAL INTERFACE)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()
  if(NOT Boost_FOUND)
    message(FATAL_ERROR "use_CGAL_Boost_support is use whereas Boost_FOUND is false.")
    return()
  endif()
  if(NOT CGAL_Boost_USE_STATIC_LIBS AND CGAL_AUTO_LINK_ENABLED)
    target_compile_definitions(${target} ${keyword} BOOST_ALL_DYN_LINK=1)
  endif()
  target_include_directories(${target} SYSTEM ${keyword} ${Boost_INCLUDE_DIRS})
  target_link_libraries(${target} ${keyword} ${Boost_LIBRARIES})
endfunction()
