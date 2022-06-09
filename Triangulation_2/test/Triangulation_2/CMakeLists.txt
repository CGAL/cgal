# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project( Triangulation_2 )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )

    include_directories (BEFORE "include")

create_single_source_cgal_program( "issue_3447.cpp" )
create_single_source_cgal_program( "issue_4025.cpp" )
create_single_source_cgal_program( "issue_4405.cpp" )
create_single_source_cgal_program( "issue_5055.cpp" )
create_single_source_cgal_program( "issue_5080.cpp" )
create_single_source_cgal_program( "test_cdt_2_projection_traits_special_case.cpp" )
create_single_source_cgal_program( "test_cdt_degenerate_case.cpp" )
create_single_source_cgal_program( "test_CDT_plus_2_with_Projection_traits_xy_3.cpp" )
create_single_source_cgal_program( "test_const_del_triangulation_2.cpp" )
create_single_source_cgal_program( "test_constrained_triangulation_2.cpp" )
create_single_source_cgal_program( "test_constrained_triangulation_with_info_2.cpp" )
create_single_source_cgal_program( "test_const_triang_plus_2.cpp" )
create_single_source_cgal_program( "test_delaunay_hierarchy_2.cpp" )
create_single_source_cgal_program( "test_delaunay_triangulation_2.cpp" )
create_single_source_cgal_program( "test_delaunay_triangulation_proj.cpp" )
create_single_source_cgal_program( "test_deprecated_projection_traits.cpp" )
create_single_source_cgal_program( "test_regular_hierarchy_2.cpp" )
create_single_source_cgal_program( "test_regular_triangulation_2.cpp" )
create_single_source_cgal_program( "test_structural_filtering_traits.cpp" )
create_single_source_cgal_program( "test_triangulation_2_bis.cpp" )
create_single_source_cgal_program( "test_triangulation_2.cpp" )
create_single_source_cgal_program( "test_triangulation_geom_traits.cpp" )
create_single_source_cgal_program( "test_triangulation_with_zip_iterator_2.cpp" )

