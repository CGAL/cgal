include(FindPackageHandleStandardArgs)
include(CheckCXXSourceCompiles)

# Usage:
#   try_CERES_with_pthread(<result_var_name> [additional linker args...])
function(try_CERES_with_pthread result_var)
    set(CERES_try_ts_source "
          #include <ceres/ceres.h>
          int main() { ceres::Problem problem; }
        ")
    set(CMAKE_REQUIRED_LIBRARIES ${ALL_CERES_LIBRARIES} ${ARGN})
    set(CMAKE_REQUIRED_INCLUDES ${CERES_INCLUDE_DIR})
    check_cxx_source_compiles("${CERES_try_ts_source}" ${result_var})
    set(${result_var} ${${result_var}} PARENT_SCOPE)
endfunction(try_CERES_with_pthread)

unset(CERES_FOUND)

find_path(CERES_INCLUDE_DIR
  NAMES
    ceres/ceres.h
  HINTS
    /usr/include/
    /usr/local/include/
)

find_library(CERES_LIBRARY
  NAMES
    ceres
  HINTS
    /usr/lib/
    /usr/lib64/
    /usr/local/lib/
    /usr/local/lib64/
)

find_package_handle_standard_args(ceres DEFAULT_MSG CERES_INCLUDE_DIR CERES_LIBRARY)

if(CERES_FOUND)
  set(CERES_INCLUDE_DIRS ${CERES_INCLUDE_DIR})
  set(CERES_LIBRARIES ${CERES_LIBRARY})

  if(UNIX AND NOT APPLE)
    try_CERES_with_pthread(CERES_without_pthread)
    if(NOT CERES_without_pthread)
      # Then check with -pthread
      try_CERES_with_pthread(CERES_with_pthread -pthread)
      if(CERES_with_pthread)
        list(APPEND CERES_LIBRARY general -pthread)
      endif(CERES_with_pthread)
    endif(NOT CERES_without_pthread)
  endif(UNIX AND NOT APPLE)

endif(CERES_FOUND)

