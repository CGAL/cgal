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
         )

find_path(LASZIP_INCLUDE_DIR
          NAMES mydefs.hpp
          PATHS /usr/local/include/LASzip/
                ${LASLIB_INCLUDE_DIR}/../LASzip/
         )

find_library(LASLIB_LIBRARIES
             NAMES las
             PATHS ENV LD_LIBRARY_PATH
                   ENV LIBRARY_PATH
                   /usr/local/lib
                   ${LASLIB_INCLUDE_DIR}/../../lib
            )

if(LASLIB_LIBRARIES)
  set(LASLIB_FOUND TRUE)
  set(LASLIB_USE_FILE "UseLASLIB")
endif()

