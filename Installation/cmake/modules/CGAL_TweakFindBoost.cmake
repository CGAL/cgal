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
  if(DEFINED CGAL_Boost_USE_STATIC_LIBS)
    # If the option is loaded from CGALConfig.h, use its value as default
    # value.  But the user will still have the choice to change the
    # value. That means that we can build the CGAL libraries using static
    # or shared Boost libraries, and after build programs using CGAL with a
    # different setting for Boost libraries.
    set(CGAL_Boost_USE_STATIC_LIBS_DEFAULT ${CGAL_Boost_USE_STATIC_LIBS})
  else()
    # Else the option is OFF by default. That means the use of shared Boost
    # libraries is the default.
    set(CGAL_Boost_USE_STATIC_LIBS_DEFAULT OFF)
  endif()

  option(CGAL_Boost_USE_STATIC_LIBS "Link with static Boost libraries" ${CGAL_Boost_USE_STATIC_LIBS_DEFAULT})
  mark_as_advanced(CGAL_Boost_USE_STATIC_LIBS)

  if(CGAL_Boost_USE_STATIC_LIBS) 
    set(Boost_USE_STATIC_LIBS ON)
  else()
    set(Boost_USE_STATIC_LIBS OFF)
    if(CGAL_AUTO_LINK_ENABLED)
      # One must add -DBOOST_ALL_DYN_LINK to DEFINITIONS to use Boost
      # auto-link with shared libraries.

      # First, add the variable to cache, if it was loaded from CGALConfig.cmake
      cache_set(CGAL_3RD_PARTY_DEFINITIONS "${CGAL_3RD_PARTY_DEFINITIONS}")
      # Then amend it
      add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS -DBOOST_ALL_DYN_LINK)
    endif()
  endif()

  set(Boost_ADDITIONAL_VERSIONS
    "1.56.1" "1.56.0" "1.56" 
    "1.55.1" "1.55.0" "1.55" 
    "1.54.1" "1.54.0" "1.54" 
    "1.53.1" "1.53.0" "1.53" 
    "1.52.1" "1.52.0" "1.52"
    "1.51.1" "1.51.0" "1.51" 
    "1.50.1" "1.50.0" "1.50" 
    "1.49.1" "1.49.0" "1.49" 
    "1.48.1" "1.48.0" "1.48" 
    "1.47.1" "1.47.0" "1.47"
    "1.46.1" "1.46.0" "1.46"
    "1.45.1" "1.45.0" "1.45"
    "1.44.1" "1.44.0" "1.44"
    "1.43.1" "1.43.0" "1.43"
    "1.42.1" "1.42.0" "1.42"
    "1.41.1" "1.41.0" "1.41"
    "1.40.1" "1.40.0" "1.40"
    "1.39.1" "1.39.0" "1.39" 
    "1.38.1" "1.38.0" "1.38" 
    "1.37.1" "1.37.0" "1.37")


  set(CGAL_TweakFindBoost)
endif()
