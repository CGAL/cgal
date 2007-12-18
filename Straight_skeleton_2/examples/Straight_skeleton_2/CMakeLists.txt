project (CGAL_STRAIGHT_SKELETON CXX)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../../cmake/modules" )

find_package(CGAL)

include(UseCGAL)

set( sources "Straight_skeleton_2.cpp" )

# Make sure the compiler can find local include files
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include")
  include_directories (BEFORE include)
endif()

# This contains an utility macro to create an executable 
# properly linked against CGAL for each .cpp in a list of files
include(MacroCreateSingleSourceCGALProgram)

create_single_source_cgal_program( ${sources} )


