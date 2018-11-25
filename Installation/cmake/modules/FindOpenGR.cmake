if ( NOT (OpenGR_LIBRARIES AND OpenGR_INCLUDE_DIR) )

# look for headers
  if (NOT OpenGR_INCLUDE_DIR)
    find_path(OpenGR_INC_DIR
              NAMES "super4pcs/shared4pcs.h")

    if(NOT OpenGR_INC_DIR)
      message(STATUS "Can not find OpenGR include directory")
    else()
      message(STATUS "Found OpenGR headers in ${OpenGR_INC_DIR}")
      set(OpenGR_INCLUDE_DIR "${OpenGR_INC_DIR}" CACHE PATH "Path to OpenGR headers files")
    endif()
  endif()

# look for library
  if (NOT OpenGR_LIBRARIES)
    find_library(OpenGR_algo_LIBRARY
                 NAMES opengr_algo)

    if (NOT OpenGR_algo_LIBRARY)
      message(STATUS "Can not find OpenGR libraries")
    else()
      message(STATUS "Found OpenGR algorithm library: ${OpenGR_algo_LIBRARY}")
      set(OpenGR_LIBRARIES "${OpenGR_algo_LIBRARY}" CACHE PATH "OpenGR libraries")
    endif()
  endif()

endif()

if ( OpenGR_LIBRARIES AND OpenGR_INCLUDE_DIR )
  set(OpenGR_FOUND ON)
endif()
