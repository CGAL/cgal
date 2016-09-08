if(CGAL_SetupCGALDependencies_included)
  return()
endif()
set(CGAL_SetupCGALDependencies_included TRUE)

if(NOT CGAL_DISABLE_GMP)
  include(CGAL_SetupGMP)
endif()

if(WITH_LEDA)
  include(CGAL_SetupLEDA)
endif()

include(CGAL_SetupBoost)

function(CGAL_setup_CGAL_dependencies target)
  if(ARGV1 STREQUAL INTERFACE)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()
  if(NOT CGAL_DISABLE_GMP)
    use_CGAL_GMP_support(${target} ${keyword})
    set(CGAL_USE_GMP TRUE CACHE INTERNAL "CGAL library is configured to use GMP")
  endif()

  if(WITH_LEDA)
    use_CGAL_LEDA_support(${target} ${keyword})
  endif()

  use_CGAL_Boost_support(${target} ${keyword})
endfunction()
