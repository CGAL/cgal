project (CGAL_STRAIGHT_SKELETON CXX)

# Allows else(), endif() etc... to ommit the expression
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

#
# Find and load CGAL settings.
#
# These lines must be kept in order
#

# 1.Guesses where FindCGAL.cmake could be
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/modules)

# 2.Try to find CGAL (not strictly required, see below)
find_package(CGAL)

# 3.Guesses where UseCGAL.cmake could be based on the found CGAL package.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CGAL_CMAKE_MODULE_PATH})

#4.Setup the parameters needed to use CGAL. UseCGAL.cmake might be found even if the CGAL package is not.
include(UseCGAL)

#
# Build the program
#

# This contains an utility macro to create an executable 
# properly linked against CGAL for each .cpp in a list of files
include(MacroCreateSingleSourceCGALProgram)

set( sources "Straight_skeleton_2.cpp" )

# Make sure the compiler can find local include files
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include")
  include_directories (BEFORE include)
endif()
 
create_single_source_cgal_program( ${sources} )


