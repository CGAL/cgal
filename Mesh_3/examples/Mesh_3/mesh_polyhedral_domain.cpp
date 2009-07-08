#include <CGAL/AABB_intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

int main()
{
  // Create polyhedron
  Polyhedron polyhedron;
  std::ifstream input("data/elephant.off");
  input >> polyhedron;

  // Create domain
  Mesh_domain domain(polyhedron);

  // Set mesh criteria
  Facet_criteria facet_criteria(25, 0.15, 0.008); // angle, size, approximation
  Cell_criteria cell_criteria(4, 0.2); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();

  // Change tetrahedron size
  Cell_criteria new_cell_criteria(4, 0.03); // radius-edge ratio, size
  Mesh_criteria new_criteria(facet_criteria, new_cell_criteria);

  // Mesh refinement
  CGAL::refine_mesh_3(c3t3, domain, new_criteria);

  // Output
  medit_file.open("out_1.mesh");
  c3t3.output_to_medit(medit_file);

  return 0;
}
