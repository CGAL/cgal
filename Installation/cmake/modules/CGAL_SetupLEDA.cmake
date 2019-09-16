#.rst:
# CGAL_SetupLEDA
# --------------
#
# The module searchs for the `LEDA` headers and library, by calling
#
# .. code-block:: cmake
#
#    find_package(LEDA)
#
# and defines the function :command:`use_CGAL_LEDA_support`.

if(CGAL_SetupLEDA_included)
  return()
endif()
set(CGAL_SetupLEDA_included TRUE)

find_package(LEDA)

#.rst:
# Provided Functions
# ^^^^^^^^^^^^^^^^^^
#
# .. command:: use_CGAL_LEDA_support
#
#    Link the target with the `LEDA` libraries::
#
#      use_CGAL_LEDA_support( target [INTERFACE] )
#
#    If the option ``INTERFACE`` is passed, the dependencies are
#    added using :command:`target_link_libraries` with the ``INTERFACE``
#    keyword, or ``PUBLIC`` otherwise.

function(use_CGAL_LEDA_support target)
  if(ARGV1 STREQUAL INTERFACE)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()
  if(NOT LEDA_FOUND)
    message(FATAL_ERROR "use_CGAL_LEDA_support is use whereas LEDA_FOUND is false.")
    return()
  endif()
  separate_arguments(LIST_LEDA_CXX_FLAGS UNIX_COMMAND "${LEDA_CXX_FLAGS}")
  separate_arguments(LIST_LEDA_DEFINITIONS UNIX_COMMAND "${LEDA_DEFINITIONS} -DCGAL_USE_LEDA")

  if(CMAKE_VERSION VERSION_LESS 3.3)
    target_compile_options(${target} ${keyword} ${LIST_LEDA_CXX_FLAGS})
  else()
    target_compile_options(${target} ${keyword} $<$<COMPILE_LANGUAGE:CXX>:${LIST_LEDA_CXX_FLAGS}>)
  endif()
  target_compile_options(${target} ${keyword} ${LIST_LEDA_DEFINITIONS})

  target_include_directories(${target} SYSTEM ${keyword} ${LEDA_INCLUDE_DIR})
  target_link_libraries(${target} ${keyword} ${LEDA_LIBRARIES} ${LEDA_LINKER_FLAGS})
endfunction()
