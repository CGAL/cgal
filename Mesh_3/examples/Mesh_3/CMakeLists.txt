# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Mesh_3_example )

cmake_minimum_required(VERSION 2.8.11)


include_directories(../../include)
include_directories(../../../Triangulation_3/include)
include_directories(../../../STL_Extension/include)
include_directories(../../../AABB_tree/include)
add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if ( MESH_3_VERBOSE )
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

find_package(CGAL COMPONENTS ImageIO)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  find_package(Boost)

  # Activate concurrency ? (turned OFF by default)
  option(ACTIVATE_CONCURRENT_MESH_3
    "Activate parallelism in Mesh_3"
    OFF)
    
  # And add -DCGAL_CONCURRENT_MESH_3 if that option is ON
  if( ACTIVATE_CONCURRENT_MESH_3 OR ENV{ACTIVATE_CONCURRENT_MESH_3} )
    add_definitions( -DCGAL_CONCURRENT_MESH_3 )
    find_package( TBB REQUIRED )
  else( ACTIVATE_CONCURRENT_MESH_3 OR ENV{ACTIVATE_CONCURRENT_MESH_3} )
    option( LINK_WITH_TBB
      "Link with TBB anyway so we can use TBB timers for profiling"
      ON)
    if( LINK_WITH_TBB )
      find_package( TBB )
    endif( LINK_WITH_TBB )
  endif()
  
  if( TBB_FOUND )
    include(${TBB_USE_FILE})
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
  endif()


  if ( Boost_FOUND AND Boost_VERSION GREATER 103400 )
    include( CGAL_CreateSingleSourceCGALProgram )

    # Compilable examples
    create_single_source_cgal_program( "mesh_implicit_sphere.cpp" )
    create_single_source_cgal_program( "mesh_implicit_sphere_variable_size.cpp" )
    create_single_source_cgal_program( "mesh_two_implicit_spheres_with_balls.cpp" )
    create_single_source_cgal_program( "mesh_implicit_domains_2.cpp" "implicit_functions.cpp" )
    create_single_source_cgal_program( "mesh_cubes_intersection.cpp" )
    create_single_source_cgal_program( "mesh_cubes_intersection_with_features.cpp" )
    create_single_source_cgal_program( "mesh_implicit_domains.cpp" "implicit_functions.cpp" )
   create_single_source_cgal_program( "mesh_polyhedral_domain.cpp" )
    create_single_source_cgal_program( "mesh_polyhedral_domain_with_features.cpp" )
    if( WITH_CGAL_ImageIO )
      if( CGAL_ImageIO_USE_ZLIB )
        create_single_source_cgal_program( "mesh_optimization_example.cpp" )
        create_single_source_cgal_program( "mesh_optimization_lloyd_example.cpp" )
        create_single_source_cgal_program( "mesh_3D_image.cpp" )
        create_single_source_cgal_program( "mesh_3D_image_variable_size.cpp" )
      else()
        message( STATUS "NOTICE: The examples mesh_3D_image.cpp, mesh_3D_image_variable_size.cpp, mesh_optimization_example.cpp and mesh_optimization_lloyd_example.cpp need CGAL_ImageIO to be configured with ZLIB support, and will not be compiled." )
      endif()
    else()
      message( STATUS "NOTICE: Some examples need the CGAL_ImageIO library, and will not be compiled." )
    endif()
#    create_single_source_cgal_program( "mesh_polyhedral_implicit_function.cpp" )
#    create_single_source_cgal_program( "mesh_polyhedral_surface_tolerance_region.cpp" )
#    create_single_source_cgal_program( "mesh_polyhedral_edge_tolerance_region.cpp" )

  else()
    message(STATUS "NOTICE: This program requires Boost >= 1.34.1, and will not be compiled.")
  endif()

else()
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

