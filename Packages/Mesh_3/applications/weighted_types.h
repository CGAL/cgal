#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// vertex and cell bases
#include <CGAL/Triangulation_vertex_base_with_info_3.h> // for Slivers_exuder.h
#include <CGAL/Complex_2_in_triangulation_vertex_base_3.h>
#include <CGAL/Complex_2_in_triangulation_cell_base_3.h>
#include <CGAL/Mesh_3/Complex_2_in_triangulation_cell_base_3.h>

// c2t3
#include <CGAL/Complex_2_in_triangulation_3.h>

// regular
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>

// traits class for reading meshes
#include <CGAL/Weighted_point_with_surface_index_geom_traits.h>

// traits class
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Regular_traits;
typedef CGAL::Weighted_point_with_surface_index_geom_traits<Regular_traits> My_traits;

// vertex and cell types
typedef CGAL::Triangulation_vertex_base_with_info_3<bool, My_traits> Vb1;
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<My_traits, Vb1> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<My_traits> Cb1;
typedef CGAL::Complex_2_in_triangulation_cell_base_3<My_traits, Cb1> Cb2;
typedef CGAL::Mesh_3::Complex_2_in_triangulation_cell_base_3<My_traits, Cb2> Cb;

// triangulation
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<My_traits, Tds> Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2T3;
