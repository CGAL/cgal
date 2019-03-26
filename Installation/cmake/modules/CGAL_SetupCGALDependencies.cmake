#.rst:
# CGAL_SetupCGALDependencies
# --------------------------
#
# The module searchs for the dependencies of the CGAL library:
#   - the `GMP/MPFR` couple,
#   - `LEDA` (optional)
#   - the `Boost` libraries (mostly the header-only libraries)
#
# and defines the variable :variable:`CGAL_FOUND` and the function
# :command:`CGAL_setup_CGAL_dependencies`.
#
# Module Input Variables
# ^^^^^^^^^^^^^^^^^^^^^^
# .. variable:: CGAL_DISABLE_GMP
#
#    If set, the `GMP` library will not be used. If
#    :variable:`WITH_LEDA` is not used either, a efficient exact
#    number types are used by CGAL kernels for exact computation.
#
# .. variable:: WITH_LEDA
#
#    If set, the `LEDA` library will be searched and used to provide
#    the exact number types used by CGAL kernels.
#
# .. variable:: CGAL_HEADER_ONLY
#
#    Set this variable if you are using the CGAL libraries as
#    header-only libraries.
#
if(CGAL_SetupCGALDependencies_included)
  return()
endif()
set(CGAL_SetupCGALDependencies_included TRUE)

#.rst:
# Used Modules
# ^^^^^^^^^^^^
#   - :module:`CGAL_SetupGMP`
if(NOT CGAL_DISABLE_GMP)
  include(${CMAKE_CURRENT_LIST_DIR}/CGAL_SetupGMP.cmake)
endif()

#.rst:
#   - :module:`CGAL_SetupLEDA`
if(WITH_LEDA)
  include(${CMAKE_CURRENT_LIST_DIR}/CGAL_SetupLEDA.cmake)
endif()

#.rst:
#   - :module:`CGAL_SetupBoost`
include(${CMAKE_CURRENT_LIST_DIR}/CGAL_SetupBoost.cmake)

#.rst:
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# .. variable:: CGAL_FOUND
#
#    Set to `TRUE` if the dependencies of CGAL were found.
if(Boost_FOUND)
  set(CGAL_FOUND TRUE)
  set_property(GLOBAL PROPERTY CGAL_FOUND TRUE)
endif()

#.rst:
#
# Provided Functions
# ^^^^^^^^^^^^^^^^^^
#
# .. command:: CGAL_setup_CGAL_dependencies
#
#   Link the target with the dependencies of CGAL::
#
#     CGAL_setup_CGAL_dependencies( target [INTERFACE] )
#
#   If the option ``INTERFACE`` is passed, the dependencies are
#   added using :command:`target_link_libraries` with the ``INTERFACE``
#   keyword, or ``PUBLIC`` otherwise.
#
function(CGAL_setup_CGAL_dependencies target)
  if(ARGV1 STREQUAL INTERFACE)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()
  if(CGAL_DISABLE_GMP)
    target_compile_definitions(${target} ${keyword} CGAL_DISABLE_GMP=1)    
  else()
    use_CGAL_GMP_support(${target} ${keyword})
    set(CGAL_USE_GMP  TRUE CACHE INTERNAL "CGAL library is configured to use GMP")
    set(CGAL_USE_MPFR TRUE CACHE INTERNAL "CGAL library is configured to use MPFR")
  endif()

  if(WITH_LEDA)
    use_CGAL_LEDA_support(${target} ${keyword})
  endif()
  
  if (CGAL_HEADER_ONLY)
    target_compile_definitions(${target} ${keyword} CGAL_HEADER_ONLY=1)
  endif()
  if (RUNNING_CGAL_AUTO_TEST)
    target_compile_definitions(${target} ${keyword} CGAL_TEST_SUITE=1)
  endif()

  use_CGAL_Boost_support(${target} ${keyword})

  foreach(dir ${CGAL_INCLUDE_DIRS})
    target_include_directories(${target} ${keyword}
      $<BUILD_INTERFACE:${dir}>)
  endforeach()
  target_include_directories(${target} ${keyword}
    $<INSTALL_INTERFACE:include>)

  # Now setup compilation flags
  if(MSVC)
    target_compile_options(${target} ${keyword}
      "-D_SCL_SECURE_NO_DEPRECATE;-D_SCL_SECURE_NO_WARNINGS"
      "/fp:strict"
      "/fp:except-"
      "/wd4503"  # Suppress warnings C4503 about "decorated name length exceeded"
      "/bigobj"  # Use /bigobj by default
      )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message( STATUS "Using Intel Compiler. Adding -fp-model strict" )
    if(WIN32)
      target_compile_options(${target} ${keyword} "/fp:strict")
    else()
      target_compile_options(${target} ${keyword} "-fp-model" "strict")
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "SunPro")
    message( STATUS "Using SunPro compiler, using STLPort 4." )
    target_compile_options(${target} ${keyword}
      "-features=extensions;-library=stlport4;-D_GNU_SOURCE")
    target_link_libraries(${target} ${keyword} "-library=stlport4")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if ( RUNNING_CGAL_AUTO_TEST )
      target_compile_options(${target} ${keyword} "-Wall")
    endif()
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 3)
      message( STATUS "Using gcc version 4 or later. Adding -frounding-math" )
      if(CMAKE_VERSION VERSION_LESS 3.3)
        target_compile_options(${target} ${keyword} "-frounding-math")
      else()
        target_compile_options(${target} ${keyword} "$<$<COMPILE_LANGUAGE:CXX>:-frounding-math>")
      endif()
    endif()
    if ( "${GCC_VERSION}" MATCHES "^4.2" )
      message( STATUS "Using gcc version 4.2. Adding -fno-strict-aliasing" )
      target_compile_options(${target} ${keyword} "-fno-strict-aliasing" )
    endif()
    if ( "${CMAKE_SYSTEM_PROCESSOR}" MATCHES "alpha" )
      message( STATUS "Using gcc on alpha. Adding -mieee -mfp-rounding-mode=d" )
      target_compile_options(${target} ${keyword} "-mieee -mfp-rounding-mode=d" )
    endif()
  endif()
endfunction()
