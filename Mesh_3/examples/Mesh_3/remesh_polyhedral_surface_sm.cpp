#define CGAL_MESH_3_VERBOSE 1
#define CGAL_MESH_3_STATS_THREADS 1
#define CGAL_CONCURRENT_COMPACT_CONTAINER_APPROXIMATE_SIZE 1

#include "draw_c3t3_surface.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/PLY_writer.h>


// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron> Mesh_domain;

using Kernel = K;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,CGAL::Parallel_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char*argv[])
{
  const char* fname = (argc>1)?argv[1]:"data/lion-head.off";
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
  // Mesh_criteria criteria(edge_size = 0.03,
  //                        facet_angle = 25, facet_size = 0.03, facet_distance = 0.003,
  //                        cell_radius_edge_ratio = 3, cell_size = 0.5);
  Mesh_criteria criteria(edge_size = 0.003,
                         facet_angle = 25, facet_size = 0.003, facet_distance = 0.0003,
                         cell_radius_edge_ratio = 3, cell_size = 0.05);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude(), manifold());

  // Output
  dump_c3t3(c3t3, "out");

  {
    using Mesh = CGAL::Surface_mesh<K::Point_3>;
    Mesh surface_mesh = draw_c3t3_surface(c3t3);

    std::ofstream out("out.ply");
    out.precision(17);
    write_ply(out, surface_mesh, "essai");
  }
}
