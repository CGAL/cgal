# This module setups the compiler for using the OpenMesh library.
# It assumes that find_package(OpenMesh) was already called.

include_directories ( ${OPENMESH_INCLUDE_DIR} )
add_definitions( -DNOMINMAX -D_USE_MATH_DEFINES  )
