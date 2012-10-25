# This module setups the compiler for using ESBTL library.
# It assumes that find_package(ESBTL) was already called.

include_directories( SYSTEM ${ESBTL_INCLUDE_DIR} )
