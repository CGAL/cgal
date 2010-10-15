#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// vertex and cell bases
#include <CGAL/Surface_mesher_vertex_base_3.h>
#include <CGAL/Surface_mesher_cell_base_3.h>
#include <CGAL/Volume_mesher_cell_base_3.h>

// c2t3
#include <CGAL/Complex_2_in_triangulation_3.h>

// delaunay
#include <CGAL/Delaunay_triangulation_3.h>

// traits class for reading meshes
#include <CGAL/Point_with_surface_index_geom_traits.h>

// traits class
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_with_surface_index_geom_traits<K> My_traits;

// vertex and cell types
typedef CGAL::Surface_mesh_vertex_base_3<My_traits> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<My_traits> Cb1;
typedef CGAL::Surface_mesh_cell_base_3<My_traits, Cb1> Cb2;
typedef CGAL::Volume_mesher_cell_base_3<My_traits, Cb2> Cb;

// triangulation
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<My_traits, Tds> Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2T3;
