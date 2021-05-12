#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// regular
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>

// IO.h must be included before vertex and cell bases.
#include <CGAL/Mesh_3/IO.h>

// vertex and cell bases
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Volume_mesher_cell_base_3.h>

// c2t3
#include <CGAL/Complex_2_in_triangulation_3.h>

// traits class for reading meshes
#include <CGAL/Weighted_point_with_surface_index_geom_traits.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

// traits class
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_filtered_traits_3<K> Regular_traits;
typedef CGAL::Weighted_point_with_surface_index_geom_traits<Regular_traits> My_traits2;
typedef CGAL::Robust_weighted_circumcenter_traits_3<My_traits2> My_traits;

// vertex and cell types
typedef CGAL::Surface_mesh_vertex_base_3<My_traits> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<My_traits> Cb1;
typedef CGAL::Surface_mesh_cell_base_3<My_traits, Cb1> Cb2;
typedef CGAL::Volume_mesher_cell_base_3<My_traits, Cb2> Cb;

// triangulation
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<My_traits, Tds> Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2T3;

const std::string format_cgal_description =
        "a file of format produced by the output operator of\n"
"  CGAL::Triangulation_3, with points\n"
"  CGAL::Weighted_point_with_surface_index and cells\n"
"  CGAL::Volume_mesher_cell_base_3.\n";
