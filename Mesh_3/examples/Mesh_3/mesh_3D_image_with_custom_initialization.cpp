#include "random_labeled_image.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Mesh_3/initialize_triangulation_from_labeled_image.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

#include <CGAL/Mesh_3/Dump_c3t3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main()
{
  /// [Create the image]
  CGAL::Image_3 image = random_labeled_image();
  /// [Create the image]

  // Domain
  Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image);

  // Mesh criteria
  Mesh_criteria criteria(facet_angle=30, facet_size=3, facet_distance=1,
                         cell_radius_edge_ratio=3, cell_size=3);

  /// [Meshing]
  C3t3 c3t3;
  initialize_triangulation_from_labeled_image(c3t3,
                                              domain,
                                              image,
                                              criteria,
                                              (unsigned char)0);
  CGAL::refine_mesh_3(c3t3, domain, criteria);
  /// [Meshing]

  // Output
  CGAL::dump_c3t3(c3t3, "out");

  return 0;
}
