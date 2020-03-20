# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Mesh_3_Examples )



add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if ( MESH_3_VERBOSE )
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

find_package(CGAL COMPONENTS ImageIO)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  find_package(Boost)

  # Activate concurrency ? (turned OFF by default)
  option(CGAL_ACTIVATE_CONCURRENT_MESH_3
    "Activate parallelism in Mesh_3"
    OFF)
    
  # And add -DCGAL_CONCURRENT_MESH_3 if that option is ON
  if( CGAL_ACTIVATE_CONCURRENT_MESH_3 OR ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3} )
    add_definitions( -DCGAL_CONCURRENT_MESH_3 )
    find_package( TBB REQUIRED )
  else( CGAL_ACTIVATE_CONCURRENT_MESH_3 OR ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3} )
    option( LINK_WITH_TBB
      "Link with TBB anyway so we can use TBB timers for profiling"
      ON)
    if( LINK_WITH_TBB )
      find_package( TBB )
    endif( LINK_WITH_TBB )
  endif()
  
  # Use Eigen
  find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
  if (NOT EIGEN3_FOUND)
    message(STATUS "This project requires the Eigen library, and will not be compiled.")
    return()
  endif()
  
  find_package(VTK QUIET COMPONENTS vtkImagingGeneral vtkIOImage NO_MODULE)
  if(VTK_FOUND)
    if(VTK_USE_FILE)
      include(${VTK_USE_FILE})
    endif()
    if ("${VTK_VERSION_MAJOR}" GREATER "5" OR VTK_VERSION VERSION_GREATER 5)
      message(STATUS "VTK found")
      if(TARGET VTK::IOImage)
        set(VTK_LIBRARIES VTK::ImagingGeneral VTK::IOImage)
      endif()
    else()
      message(STATUS "VTK version 6.0 or greater is required")
      endif()
    else()
      message(STATUS "VTK was not found")
    endif()


    # Compilable examples
    create_single_source_cgal_program( "mesh_hybrid_mesh_domain.cpp" )
    CGAL_target_use_Eigen(mesh_hybrid_mesh_domain)
    
    create_single_source_cgal_program( "mesh_implicit_sphere.cpp" )
    CGAL_target_use_Eigen(mesh_implicit_sphere)
    
    create_single_source_cgal_program( "mesh_implicit_sphere_variable_size.cpp" )
    CGAL_target_use_Eigen(mesh_implicit_sphere_variable_size)
    
    create_single_source_cgal_program( "mesh_two_implicit_spheres_with_balls.cpp" )
    CGAL_target_use_Eigen(mesh_two_implicit_spheres_with_balls)
    
    create_single_source_cgal_program( "mesh_implicit_domains_2.cpp" "implicit_functions.cpp" )
    CGAL_target_use_Eigen(mesh_implicit_domains_2)
    
    create_single_source_cgal_program( "mesh_cubes_intersection.cpp" )
    CGAL_target_use_Eigen(mesh_cubes_intersection)
    
    create_single_source_cgal_program( "mesh_cubes_intersection_with_features.cpp" )
    CGAL_target_use_Eigen(mesh_cubes_intersection_with_features)
    
    create_single_source_cgal_program( "mesh_implicit_domains.cpp" "implicit_functions.cpp" )
    CGAL_target_use_Eigen(mesh_implicit_domains)
    
    create_single_source_cgal_program( "mesh_polyhedral_domain.cpp" )
    CGAL_target_use_Eigen(mesh_polyhedral_domain)
    
    create_single_source_cgal_program( "mesh_polyhedral_domain_sm.cpp" )
    CGAL_target_use_Eigen(mesh_polyhedral_domain_sm)
    
    create_single_source_cgal_program( "mesh_polyhedral_domain_with_surface_inside.cpp" )
    CGAL_target_use_Eigen(mesh_polyhedral_domain_with_surface_inside)
    
    create_single_source_cgal_program( "remesh_polyhedral_surface.cpp" )
    CGAL_target_use_Eigen(remesh_polyhedral_surface)
    
    create_single_source_cgal_program( "remesh_polyhedral_surface_sm.cpp" )
    CGAL_target_use_Eigen(remesh_polyhedral_surface_sm)
    
    create_single_source_cgal_program( "mesh_polyhedral_domain_with_features.cpp" )
    CGAL_target_use_Eigen(mesh_polyhedral_domain_with_features)
    
    create_single_source_cgal_program( "mesh_polyhedral_domain_with_features_sm.cpp" )
    CGAL_target_use_Eigen(mesh_polyhedral_domain_with_features_sm)
    
    create_single_source_cgal_program( "mesh_polyhedral_domain_with_lipschitz_sizing.cpp" )
    CGAL_target_use_Eigen(mesh_polyhedral_domain_with_lipschitz_sizing)
    
    create_single_source_cgal_program( "mesh_polyhedral_complex.cpp" )
    CGAL_target_use_Eigen(mesh_polyhedral_complex)
    
    create_single_source_cgal_program( "mesh_polyhedral_complex_sm.cpp" )
    CGAL_target_use_Eigen(mesh_polyhedral_complex_sm)
    
    if( WITH_CGAL_ImageIO )
      if( VTK_FOUND AND ("${VTK_VERSION_MAJOR}" GREATER "5" OR VTK_VERSION VERSION_GREATER 5) )
        add_executable ( mesh_3D_gray_vtk_image mesh_3D_gray_vtk_image.cpp )
        CGAL_target_use_Eigen(mesh_3D_gray_vtk_image)
        target_link_libraries(  mesh_3D_gray_vtk_image ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${VTK_LIBRARIES})
        cgal_add_test( mesh_3D_gray_vtk_image )
        add_to_cached_list( CGAL_EXECUTABLE_TARGETS mesh_3D_gray_vtk_image )
      endif()
 
      create_single_source_cgal_program( "mesh_3D_gray_image.cpp" )
      CGAL_target_use_Eigen(mesh_3D_gray_image)
      
      create_single_source_cgal_program( "mesh_3D_gray_image_multiple_values.cpp" )
      CGAL_target_use_Eigen(mesh_3D_gray_image_multiple_values)
      
      create_single_source_cgal_program( "mesh_3D_image_with_features.cpp" )
      CGAL_target_use_Eigen(mesh_3D_image_with_features)
      
      if( CGAL_ImageIO_USE_ZLIB )
        create_single_source_cgal_program( "mesh_optimization_example.cpp" )
        CGAL_target_use_Eigen(mesh_optimization_example)

        create_single_source_cgal_program( "mesh_optimization_lloyd_example.cpp" )
        CGAL_target_use_Eigen(mesh_optimization_lloyd_example)

        create_single_source_cgal_program( "mesh_3D_image.cpp" )
        CGAL_target_use_Eigen(mesh_3D_image)
        
        create_single_source_cgal_program( "mesh_3D_image_with_custom_initialization.cpp" )
        CGAL_target_use_Eigen(mesh_3D_image_with_custom_initialization)

        create_single_source_cgal_program( "mesh_3D_image_variable_size.cpp" )
        CGAL_target_use_Eigen(mesh_3D_image_variable_size)
      else()
        message( STATUS "NOTICE: The examples mesh_3D_image.cpp, mesh_3D_image_variable_size.cpp, mesh_optimization_example.cpp and mesh_optimization_lloyd_example.cpp need CGAL_ImageIO to be configured with ZLIB support, and will not be compiled." )
      endif()
    else()
      message( STATUS "NOTICE: Some examples need the CGAL_ImageIO library, and will not be compiled." )
    endif()
    
#    create_single_source_cgal_program( "mesh_polyhedral_implicit_function.cpp" )
#    create_single_source_cgal_program( "mesh_polyhedral_surface_tolerance_region.cpp" )
#    create_single_source_cgal_program( "mesh_polyhedral_edge_tolerance_region.cpp" )

    if(CGAL_ACTIVATE_CONCURRENT_MESH_3 AND TBB_FOUND AND TARGET ${target})
      foreach(target
        mesh_3D_image_with_features
        mesh_3D_image
        mesh_polyhedral_domain
        mesh_3D_image_with_custom_initialization
        mesh_two_implicit_spheres_with_balls
        mesh_optimization_lloyd_example
        mesh_optimization_example
        mesh_implicit_sphere
        mesh_polyhedral_complex_sm
        mesh_implicit_sphere_variable_size
        mesh_polyhedral_domain_sm
        mesh_polyhedral_domain_with_lipschitz_sizing
        mesh_polyhedral_complex
        mesh_polyhedral_domain_with_features
        mesh_3D_image_variable_size)
        CGAL_target_use_TBB(${target})
      endforeach()
    endif()

else()
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

