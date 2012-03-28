#define CGAL_MESH_3_VERBOSE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Gray_image_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <functional>

typedef short Image_word_type;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Gray_image_mesh_domain_3<CGAL::Image_3,K, 
                                       Image_word_type,
                                       std::binder1st< std::less<Image_word_type> > > Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char** argv)
{
  // Loads image
  CGAL::Image_3 image;
  if(!image.read(argv[1])) return 1;

  // Domain
  Mesh_domain domain(image, std::bind1st(std::less<Image_word_type>(), -548), -10000);

  // Mesh criteria
  Mesh_criteria criteria(facet_angle=30, facet_size=1, facet_distance=0.05,
                         cell_radius_edge_ratio=100, cell_size=0);

  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, 
                         CGAL::parameters::no_exude(), CGAL::parameters::no_perturb());

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);

  return 0;
}
