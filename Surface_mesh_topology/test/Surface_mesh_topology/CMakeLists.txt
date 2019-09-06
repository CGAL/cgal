project(Surface_mesh_topology_Tests)

cmake_minimum_required(VERSION 3.1)

find_package(CGAL COMPONENTS Qt5)
if(CGAL_Qt5_FOUND)
    add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

if ( CGAL_FOUND )
  create_single_source_cgal_program("path_tests.cpp")
  create_single_source_cgal_program("tests_path_on_surface.cpp")
  create_single_source_cgal_program("test_homotopy.cpp")
  create_single_source_cgal_program("fundamental_group_of_the_circle.cpp")
  create_single_source_cgal_program("fundamental_group_of_the_torus.cpp")
  create_single_source_cgal_program("homotopy_big_cylinder.cpp")
  create_single_source_cgal_program("homotopy_double_torus_with_holes.cpp")
  create_single_source_cgal_program("homotopy_rond_point_saucisse.cpp")
  create_single_source_cgal_program("homotopy_small_cylinder.cpp")
  create_single_source_cgal_program("path_with_rle_deformation_tests.cpp")

  if(CGAL_Qt5_FOUND)
    target_link_libraries(path_with_rle_deformation_tests PUBLIC CGAL::CGAL_Qt5)
    target_link_libraries(test_homotopy PUBLIC CGAL::CGAL_Qt5)
    target_link_libraries(fundamental_group_of_the_circle PUBLIC CGAL::CGAL_Qt5)
    target_link_libraries(fundamental_group_of_the_torus PUBLIC CGAL::CGAL_Qt5)
  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()
