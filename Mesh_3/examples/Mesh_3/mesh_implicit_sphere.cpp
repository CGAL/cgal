#define CGAL_MESH_3_VERBOSE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_cell_base_3.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Delaunay_triangulation_3.h>

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point = K::Point_3;
using Function = FT(const Point&);
using Mesh_domain = CGAL::Labeled_mesh_domain_3<K>;

using Concurrency_tag =
#ifdef CGAL_CONCURRENT_MESH_3
        CGAL::Parallel_tag;
#else
        CGAL::Sequential_tag;
#endif

// Triangulation
using Dt_vb = CGAL::Triangulation_vertex_base_3<K>;
using Vb = CGAL::Mesh_vertex_base_3<K, Mesh_domain, Dt_vb>;
using Dt_cb = CGAL::Delaunay_triangulation_cell_base_3<K>;
using Cb = CGAL::Mesh_cell_base_3<K, Mesh_domain, Dt_cb>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb, Concurrency_tag>;
using Tr = CGAL::Delaunay_triangulation_3<K, Tds>;

using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

namespace params = CGAL::parameters;

// Function
FT sphere_function (const Point& p)
{ return CGAL::squared_distance(p, Point(CGAL::ORIGIN))-1; }

int main()
{
  /// [Domain creation] (Warning: Sphere_3 constructor uses squared radius !)
  Mesh_domain domain =
    Mesh_domain::create_implicit_mesh_domain( sphere_function,
                                              K::Sphere_3(CGAL::ORIGIN, K::FT(2)));
  /// [Domain creation]

  // Mesh criteria
  Mesh_criteria criteria(params::facet_angle(30).facet_size(0.1).facet_distance(0.025).
                         cell_radius_edge_ratio(2).cell_size(0.1));

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  CGAL::IO::write_MEDIT(medit_file, c3t3);
  medit_file.close();

  return 0;
}
