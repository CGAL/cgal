#define CGAL_MESH_3_STATS_THREADS 1

#include "draw_c3t3_surface.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Surface_mesh/IO/PLY.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,CGAL::Parallel_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

namespace params = CGAL::parameters;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("meshes/lion-head.off");
  // Load a polyhedron
  Polyhedron poly;
  std::ifstream input(fname);
  input >> poly;

  if (!CGAL::is_triangle_mesh(poly)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // Create a vector with only one element: the pointer to the polyhedron.
  std::vector<Polyhedron*> poly_ptrs_vector(1, &poly);

  // Create a polyhedral domain, with only one polyhedron,
  // and no "bounding polyhedron", so the volumetric part of the domain will be
  // empty.
  Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());

  // Get sharp features
  domain.detect_features(); //includes detection of borders

  // Mesh criteria
  Mesh_criteria criteria(params::edge_size(0.025).
                                 facet_angle(25).
                                 facet_size(0.1).
                                 facet_distance(0.001));

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, params::no_perturb().no_exude());

  // Output
  dump_c3t3(c3t3, "out");

#if CGAL_USE_BASIC_VIEWER
  {
    using Mesh = CGAL::Surface_mesh<K::Point_3>;
    Mesh surface_mesh = draw_c3t3_surface(c3t3);
    std::ofstream out("out.ply");
    out.precision(17);
    CGAL::IO::write_PLY(out, surface_mesh, "essai");
  }
#endif
}
