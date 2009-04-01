# This is the CMake script for compiling a CGAL application.

project( pdb ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS PDB )
find_package(ImageMagick COMPONENTS Magick++)
include( ${CGAL_USE_FILE} )

find_package ( Boost COMPONENTS program_options )


link_directories (${PACKAGE_ROOT}/../CGALPDB/src/CGALPDB)


if (CGAL_PDB_FOUND  AND  Boost_PROGRAM_OPTIONS_FOUND )

  link_libraries(${Boost_LIBRARIES} ${CGAL_PDB_LIBRARY})

  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( pdb_cat.cpp )
  create_single_source_cgal_program( pdb_echo.cpp )
  create_single_source_cgal_program( pdb_split.cpp )
  create_single_source_cgal_program( pdb_transform.cpp )
  create_single_source_cgal_program( pdb_tree.cpp )
  create_single_source_cgal_program( pdb_info.cpp )

  if (ImageMagick_Magick++_FOUND)
      include_directories(${ImageMagick_INCLUDE_DIRS})
      link_libraries(${ImageMagick_LIBRARIES})
      message(STATUS "ImageMagic headers found at "${ImageMagick_INCLUDE_DIRS})
      create_single_source_cgal_program( pdb_distance_matrix.cpp )
  else()

     message(STATUS "NOTICE: The distance matrix demo requires ImageMagick, and will not be compiled.")
  endif()

else()

  message(STATUS "NOTICE: This demo requires PDB and Boost program_options, and will not be compiled.")

endif()
