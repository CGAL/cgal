# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Mesh_3_test )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core ImageIO)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)
  include_directories (BEFORE ../../../AABB_tree/include)

  create_single_source_cgal_program( "test_backward_compatibility_MeshFoobarCriteria_3.cpp" )
  create_single_source_cgal_program( "test_boost_has_xxx.cpp" )
  create_single_source_cgal_program( "test_c3t3.cpp" )
  create_single_source_cgal_program( "test_c3t3_with_features.cpp" )
  create_single_source_cgal_program( "test_criteria.cpp" )
  create_single_source_cgal_program( "test_domain_with_polyline_features.cpp" )
  create_single_source_cgal_program( "test_mesh_criteria_creation.cpp" )
  if(CGAL_ImageIO_USE_ZLIB)
    create_single_source_cgal_program( "test_meshing_3D_image.cpp" )
  else()
    message(STATUS "test_meshing_3D_image requires the ZLIB library, and will not be compiled.")
  endif()
  create_single_source_cgal_program( "test_meshing_implicit_function.cpp" )
  create_single_source_cgal_program( "test_meshing_polyhedron.cpp" )
  create_single_source_cgal_program( "test_meshing_polyhedron_with_features.cpp" )
  create_single_source_cgal_program( "test_robust_weighted_circumcenter.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

