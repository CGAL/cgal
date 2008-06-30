project( Alpha_shapes_3 ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED)


include( ${CGAL_USE_FILE} )
include( CreateSingleSourceCGALProgram )

create_single_source_cgal_program( alpha_shapes_3.cpp )
create_single_source_cgal_program( weighted_alpha_shapes_3.cpp )

