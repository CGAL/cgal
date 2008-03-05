# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (STL_Extension_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( test_Cache.cpp )
  create_single_source_cgal_program( test_Compact_container.cpp )
  create_single_source_cgal_program( test_composition.cpp )
  create_single_source_cgal_program( test_Concatenate_iterator.cpp )
  create_single_source_cgal_program( test_Flattening_iterator.cpp )
  create_single_source_cgal_program( test_Handle_with_policy.cpp )
  create_single_source_cgal_program( test_In_place_list.cpp )
  create_single_source_cgal_program( test_lexcompare_outputrange.cpp )
  create_single_source_cgal_program( test_multiset.cpp )
  create_single_source_cgal_program( test_Nested_iterator.cpp )
  create_single_source_cgal_program( test_new_binders1.cpp )
  create_single_source_cgal_program( test_stl_extension.cpp )
  create_single_source_cgal_program( test_type_traits.cpp )
  create_single_source_cgal_program( test_vector.cpp )

endif()
