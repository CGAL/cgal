# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_modeling_ )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    include( CGAL_CreateSingleSourceCGALProgram )

    include_directories (BEFORE "../../include")

    create_single_source_cgal_program( "all_roi_assign_example.cpp" )
    create_single_source_cgal_program( "all_roi_assign_example_custom_polyhedron.cpp" )
    create_single_source_cgal_program( "all_roi_assign_example_Surface_mesh.cpp" )
    create_single_source_cgal_program( "custom_weight_for_edges_example.cpp" )
    create_single_source_cgal_program( "deform_polyhedron_with_custom_pmap_example.cpp" )
    create_single_source_cgal_program( "k_ring_roi_translate_rotate_example.cpp" )
    create_single_source_cgal_program( "k_ring_roi_translate_rotate_Surface_mesh.cpp" )

    find_package( OpenMesh QUIET )
    if ( OpenMesh_FOUND )
    include( UseOpenMesh )
      create_single_source_cgal_program( "all_roi_assign_example_with_OpenMesh.cpp" )
      target_link_libraries( all_roi_assign_example_with_OpenMesh ${OPENMESH_LIBRARIES} )
    else()
      message(STATUS "Example that use OpenMesh will not be compiled.")
    endif()

  else()
    message(STATUS "NOTICE: These examples require the Eigen library, version 3.2 or later and will not be compiled.")
  endif()
else()
  
    message(STATUS "NOTICE: These exmaples require the CGAL library, and will not be compiled.")
  
endif()

