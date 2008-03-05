# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Distance_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( tst20.cpp )
  create_single_source_cgal_program( tst21.cpp )
  create_single_source_cgal_program( tst22.cpp )
  create_single_source_cgal_program( tst23.cpp )
  create_single_source_cgal_program( tst24.cpp )
  create_single_source_cgal_program( tst25.cpp )
  create_single_source_cgal_program( tst26.cpp )
  create_single_source_cgal_program( tst27.cpp )
  create_single_source_cgal_program( tst28.cpp )
  create_single_source_cgal_program( tst29.cpp )
  create_single_source_cgal_program( tst30.cpp )
  create_single_source_cgal_program( tst31.cpp )
  create_single_source_cgal_program( tst32.cpp )
  create_single_source_cgal_program( tst33.cpp )
  create_single_source_cgal_program( tst34.cpp )

endif()
