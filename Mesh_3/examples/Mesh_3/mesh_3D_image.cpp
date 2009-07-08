#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Image_3 Image;
typedef CGAL::Labeled_image_mesh_domain_3<Image,K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;

int main()
{
  // Loads image
  Image image;
  image.read("data/liver.inr.gz");

  // Domain
  Mesh_domain domain(image);

  // Mesh criteria
  Facet_criteria facet_criteria(30, 6, 4); // angle, size, approximation
  Cell_criteria cell_criteria(3, 8); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);

  return 0;
}
