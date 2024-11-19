# - Try to find LASLIB
# Once done this will define
#
#  LASLIB_FOUND = LASLIB_FOUND - TRUE
#  LASLIB_INCLUDE_DIR - include directory for LASlib
#  LASZIP_INCLUDE_DIR - include directory for LASzip
#  LASLIB_LIBRARIES   - the libraries (as targets)

# first look in user defined locations
find_path(LASLIB_INCLUDE_DIR
          NAMES lasreader.hpp
         PATHS /usr/local/include/LASlib/
                   ENV LASLIB_INC_DIR
         )

find_path(LASZIP_INCLUDE_DIR
          NAMES mydefs.hpp
          PATHS /usr/local/include/LASzip/
                ${LASLIB_INCLUDE_DIR}/../../LASzip/src
                ${LASLIB_INCLUDE_DIR}/../LASzip
                ${LASLIB_INCLUDE_DIR}
         )

find_library(LASLIB_LIBRARIES
             NAMES las
             PATHS ENV LD_LIBRARY_PATH
                   ENV LIBRARY_PATH
                   /usr/local/lib
                   /usr/local/lib/LASlib
                   ${LASLIB_INCLUDE_DIR}/../../lib
                  ENV LASLIB_LIB_DIR
            )
if (NOT LASLIB_LIBRARIES)
  #library was renamed in recent versions of LAStools
  find_library(LASLIB_LIBRARIES
               NAMES LASlib
               PATHS ENV LD_LIBRARY_PATH
                     ENV LIBRARY_PATH
                     /usr/local/lib
                     /usr/local/lib/LASlib
                     ${LASLIB_INCLUDE_DIR}/../../lib
                    ENV LASLIB_LIB_DIR
              )
endif()

if(LASLIB_LIBRARIES AND LASLIB_INCLUDE_DIR AND LASZIP_INCLUDE_DIR)
  if (NOT ${LASLIB_INCLUDE_DIR} STREQUAL ${LASZIP_INCLUDE_DIR})
    list(APPEND LASLIB_INCLUDE_DIR ${LASZIP_INCLUDE_DIR})
  endif()
  set(LASLIB_FOUND TRUE)
  set(LASLIB_USE_FILE "UseLASLIB")
endif()
