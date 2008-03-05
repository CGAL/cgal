# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Polyhedron_IO_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( iv2off.cpp )
  create_single_source_cgal_program( off2iv.cpp )
  create_single_source_cgal_program( off2off.cpp )
  create_single_source_cgal_program( off2stl.cpp )
  create_single_source_cgal_program( off2vrml.cpp )
  create_single_source_cgal_program( off2wav.cpp )
  create_single_source_cgal_program( off_bbox.cpp )
  create_single_source_cgal_program( off_glue.cpp )
  create_single_source_cgal_program( off_transform.cpp )
  create_single_source_cgal_program( polyhedron2vrml.cpp )
  create_single_source_cgal_program( polyhedron_copy.cpp )
  create_single_source_cgal_program( terr_trian.cpp )

endif()
