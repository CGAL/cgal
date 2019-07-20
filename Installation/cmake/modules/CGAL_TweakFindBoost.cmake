# - Defines Boost_USE_STATIC_LIBS and Boost_ADDITIONAL_VERSIONS
# 
# This module sets the CMake variables:
# 
# == Boost_USE_STATIC_LIBS ==
# 
# The option CGAL_Boost_USE_STATIC_LIBS is created in the cache, as
# advanced option. If CGALConfig.cmake has been loaded, the default value
# of that option is the value loaded from CGALConfig.cmake (this file was
# created during the configuration of CGAL libraries). Otherwise, the
# default value of that option is OFF.
# 
# The variable Boost_USE_STATIC_LIBS is set to the value of the option
# CGAL_Boost_USE_STATIC_LIBS.
# 
# Additionally, if Boost_USE_STATIC_LIBS is OFF, and the auto-linking is
# enabled, the definition BOOST_ALL_DYN_LINK is added to
# CGAL_3RD_PARTY_DEFINITIONS, so that the auto-linking feature on Windows
# knows that it must search for dynamic libraries.
# 
# == Boost_ADDITIONAL_VERSIONS ==
# 
# The variable Boost_ADDITIONAL_VERSIONS is filled with a long list of
# Boost versions. That allows the module FindBoost to find more recent
# Boost versions, even if the file FindBoost.cmake is old.

if( NOT CGAL_TweakFindBoost )
  if(POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
  endif()
  if(DEFINED CGAL_Boost_USE_STATIC_LIBS)
    # If the option is loaded from CGALConfig.cmake, use its value as default
    # value.  But the user will still have the choice to change the
    # value. That means that we can build the CGAL libraries using static
    # or shared Boost libraries, and after build programs using CGAL with a
    # different setting for Boost libraries.
    set(CGAL_Boost_USE_STATIC_LIBS_DEFAULT ${CGAL_Boost_USE_STATIC_LIBS})
  else()
    # Else the option default is related to BUILD_SHARED_LIBS.
    if(BUILD_SHARED_LIBS OR NOT DEFINED BUILD_SHARED_LIBS)
      set(CGAL_Boost_USE_STATIC_LIBS_DEFAULT OFF)
    else()
      set(CGAL_Boost_USE_STATIC_LIBS_DEFAULT ON)
    endif()
  endif()

  option(Boost_DEBUG "Activate the debug messages of the script FindBoost" OFF)
  mark_as_advanced(Boost_DEBUG)

  option(CGAL_Boost_USE_STATIC_LIBS "Link with static Boost libraries" ${CGAL_Boost_USE_STATIC_LIBS_DEFAULT})
  mark_as_advanced(CGAL_Boost_USE_STATIC_LIBS)

  if(CGAL_Boost_USE_STATIC_LIBS) 
    set(Boost_USE_STATIC_LIBS ON)
  else()
    set(Boost_USE_STATIC_LIBS OFF)
    if(CGAL_AUTO_LINK_ENABLED)
      # One must add -DBOOST_ALL_DYN_LINK to DEFINITIONS to use Boost
      # auto-link with shared libraries.

      list(APPEND CGAL_3RD_PARTY_DEFINITIONS -DBOOST_ALL_DYN_LINK)
      set(CGAL_3RD_PARTY_DEFINITIONS "${CGAL_3RD_PARTY_DEFINITIONS}"
	CACHE INTERNAL "3rd party definitions for CGAL")
    endif()
  endif()

  set(Boost_ADDITIONAL_VERSIONS
    "1.89.1" "1.89.0" "1.89"
    "1.88.1" "1.88.0" "1.88"
    "1.87.1" "1.87.0" "1.87"
    "1.86.1" "1.86.0" "1.86"
    "1.85.1" "1.85.0" "1.85"
    "1.84.1" "1.84.0" "1.84"
    "1.83.1" "1.83.0" "1.83"
    "1.82.1" "1.82.0" "1.82"
    "1.81.1" "1.81.0" "1.81"
    "1.80.1" "1.80.0" "1.80"
    "1.79.1" "1.79.0" "1.79"
    "1.78.1" "1.78.0" "1.78"
    "1.77.1" "1.77.0" "1.77"
    "1.76.1" "1.76.0" "1.76"
    "1.75.1" "1.75.0" "1.75"
    "1.74.1" "1.74.0" "1.74"
    "1.73.1" "1.73.0" "1.73"
    "1.72.1" "1.72.0" "1.72"
    "1.71.1" "1.71.0" "1.71"
    "1.70.1" "1.70.0" "1.70"
    "1.69.1" "1.69.0" "1.69"
    "1.68.1" "1.68.0" "1.68"
    "1.67.1" "1.67.0" "1.67"
    "1.66.1" "1.66.0" "1.66"
    "1.65.1" "1.65.0" "1.65"
    "1.64.1" "1.64.0" "1.64"
    "1.63.1" "1.63.0" "1.63"
    "1.62.1" "1.62.0" "1.62"
    "1.61.1" "1.61.0" "1.61"
    "1.60.1" "1.60.0" "1.60"
    "1.59.1" "1.59.0" "1.59"
    "1.58.1" "1.58.0" "1.58"
    "1.57.1" "1.57.0" "1.57"
    )


  set(CGAL_TweakFindBoost ON)
endif()
