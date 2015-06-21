#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_image_mesh_domain_3<CGAL::Image_3,K> Mesh_domain;

// Triangulation
#ifdef CGAL_CONCURRENT_MESH_3
  typedef CGAL::Mesh_triangulation_3<
    Mesh_domain,
    CGAL::Kernel_traits<Mesh_domain>::Kernel, // Same as sequential
    CGAL::Parallel_tag                        // Tag to activate parallelism
  >::type Tr;
#else
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main()
{
  // Domain
  CGAL::Image_3 image;
  image.read("data/liver.inr.gz");
  Mesh_domain domain(image);

  // Mesh criteria
  Mesh_criteria criteria(facet_angle=30, facet_distance=1.2,
                         cell_radius_edge_ratio=2);

  // Mesh generation and optimization in one call
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      lloyd(time_limit=30),
                                      no_perturb(),
                                      exude(time_limit=10, sliver_bound=10));
  
  // Mesh generation and optimization in several call
  C3t3 c3t3_bis = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                          no_perturb(), no_exude());
  
  CGAL::lloyd_optimize_mesh_3(c3t3_bis, domain, time_limit=30);
  CGAL::exude_mesh_3(c3t3_bis, sliver_bound=10, time_limit=10);
  
  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);
  
  std::ofstream medit_file_bis("out_bis.mesh");
  c3t3_bis.output_to_medit(medit_file_bis);

  return 0;
}
