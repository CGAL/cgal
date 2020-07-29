#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

// Polyhedron type
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main()
{
  // Load a polyhedron
  Polyhedron poly;
  std::ifstream input("data/lion-head.off");
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
  Mesh_criteria criteria(edge_size = 0.025,
                         facet_angle = 25,
                         facet_size = 0.1,
                         facet_distance = 0.001);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  // Output the facets of the c3t3 to an OFF file. The facets will not be
  // oriented.
  std::ofstream off_file("out.off");
  c3t3.output_boundary_to_off(off_file);

  return off_file.fail() ? EXIT_FAILURE : EXIT_SUCCESS;
}
