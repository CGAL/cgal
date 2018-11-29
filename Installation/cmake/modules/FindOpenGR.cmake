if ( OpenGR_INCLUDE_DIR )

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

endif()

if ( OpenGR_INCLUDE_DIR )
  set(OpenGR_FOUND ON)
endif()
