#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

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

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

namespace params = CGAL::parameters;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("images/liver.inr.gz");
  // Domain
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image);

  // Mesh criteria
  Mesh_criteria criteria(params::facet_angle(30).facet_distance(1.2).
                                 cell_radius_edge_ratio(2));

  // Mesh generation and optimization in one call
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      params::lloyd(params::time_limit(30)).
                                      no_perturb().
                                      exude(params::time_limit(10).sliver_bound(10)));

  // Mesh generation and optimization in several call
  C3t3 c3t3_bis = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                          params::no_perturb().no_exude());

  CGAL::lloyd_optimize_mesh_3(c3t3_bis, domain, params::time_limit(30));
  CGAL::exude_mesh_3(c3t3_bis, params::sliver_bound(10), params::time_limit(10));

  // Output
  std::ofstream medit_file("out.mesh");
  CGAL::IO::write_MEDIT(medit_file, c3t3);
  medit_file.close();

  std::ofstream medit_file_bis("out_bis.mesh");
  CGAL::IO::write_MEDIT(medit_file_bis, c3t3_bis);
  medit_file_bis.close();

  return 0;
}
