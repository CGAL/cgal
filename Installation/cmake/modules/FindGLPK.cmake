# This file sets up GLPK for CMake. Once done this will define
#  GLPK_FOUND             - system has GLPK lib
#  GLPK_INCLUDE_DIR       - the GLPK include directory
#  GLPK_LIBRARIES         - Link these to use GLPK


# Is it already configured?
if (NOT GLPK_FOUND)

    # first look in user defined locations
    find_path(GLPK_INCLUDE_DIR
                NAMES glpk.h
                PATHS /usr/local/include/LASlib/
                ENV GLPK_INC_DIR
             )

    find_library(GLPK_LIBRARIES
                 NAMES libglpk glpk
                 PATHS ENV LD_LIBRARY_PATH
                       ENV LIBRARY_PATH
                       /usr/local/lib
                       ${GLPK_INCLUDE_DIR}/../lib
                      ENV GLPK_LIB_DIR
                )

    if(GLPK_LIBRARIES AND GLPK_INCLUDE_DIR)
      set(GLPK_FOUND TRUE)
    endif()

endif()
