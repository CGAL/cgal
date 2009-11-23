#ifndef C2T3_TYPE_H
#define C2T3_TYPE_H

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Mesh_3/Robust_weighted_circumcenter_filtered_traits_3.h>

// traits class
typedef CGAL::Exact_predicates_inexact_constructions_kernel K2;
typedef CGAL::Robust_weighted_circumcenter_filtered_traits_3<K2>  Traits;

// typedef CGAL::Regular_triangulation_euclidean_traits_3<K>  Traits;

// vertex and cell types
typedef CGAL::Surface_mesh_vertex_base_3<Traits> Vb;
typedef CGAL::Surface_mesh_cell_base_3<Traits> Cb;
typedef CGAL::Triangulation_cell_base_with_circumcenter_3<Traits, Cb> Cb_with_circumcenter;

// triangulation
typedef CGAL::Triangulation_data_structure_3<Vb, Cb_with_circumcenter> Tds;
typedef CGAL::Regular_triangulation_3<Traits, Tds> Tr;

typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

#endif
