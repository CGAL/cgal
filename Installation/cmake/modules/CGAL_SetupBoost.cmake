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

find_package( Boost 1.48 REQUIRED )

if(Boost_FOUND AND Boost_VERSION VERSION_LESS 1.70)
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
  if(TARGET Boost::boost)
    target_link_libraries(${target} ${keyword} Boost::boost)
  else()
    target_include_directories(${target} SYSTEM ${keyword} ${Boost_INCLUDE_DIRS})
    target_link_libraries(${target} ${keyword} ${Boost_LIBRARIES})
  endif()
endfunction()
