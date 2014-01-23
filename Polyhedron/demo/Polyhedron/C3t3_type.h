#ifndef C3T3_TYPE_H
#define C3T3_TYPE_H

#include "Polyhedron_type.h"

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Triangle_accessor_3.h>

#include <vector>

namespace c3t3_type_h {

// @TODO: Is that the right kernel?!
typedef CGAL::Polyhedral_mesh_domain_3<
  Polyhedron,
  Kernel,
  CGAL::Triangle_accessor_3<Polyhedron,Kernel>,
  CGAL::Tag_true,
  CGAL::Tag_false // Kernel has already been patched with exact intersection
  >                           Mesh_domain_without_features;
typedef CGAL::Mesh_domain_with_polyline_features_3<
  Mesh_domain_without_features>             Mesh_domain;
typedef Mesh_domain                         MD;

typedef CGAL::Robust_weighted_circumcenter_filtered_traits_3<Kernel> Geom_traits;

// Triangulation
typedef CGAL::Mesh_vertex_base_3<Geom_traits, MD>                   Vertex_base;

typedef CGAL::Compact_mesh_cell_base_3<Geom_traits, MD>             Cell_base;
typedef CGAL::Triangulation_data_structure_3<Vertex_base,Cell_base> Tds;
typedef CGAL::Regular_triangulation_3<Geom_traits, Tds>             Tr;

} // end namespace c3t3_type_h

using c3t3_type_h::Mesh_domain;

// 3D complex
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  c3t3_type_h::Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_segment_index> C3t3;

#endif // C3T3_TYPE_H
