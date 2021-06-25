#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/generate_weights_from_labeled_image.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/IO/File_binary_mesh_3.h>

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

using Word_type = unsigned char;
using Subdomain_index = int;

int main(int argc, char* argv[])
{
  /// [Loads image]
  const char* fname = (argc > 1) ? argv[1] : "data/liver.inr.gz";
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  /// [Loads image]

  /// [Generate weights]
  const float sigma = (argc > 2) ? atof(argv[2]) : 1.f;
  CGAL::Image_3 img_weights =
    CGAL::Mesh_3::generate_weights(image, sigma, (unsigned char)(1));
  /// [Generate weights]

  /// [Domain creation]
  Mesh_domain domain
    = Mesh_domain::create_labeled_image_mesh_domain(image,
                                                    weights = img_weights,
                                                    relative_error_bound = 1e-6);
  /// [Domain creation]

  // Mesh criteria
  Mesh_criteria criteria(facet_angle=30, facet_size=6, facet_distance=0.5,
                         cell_radius_edge_ratio=3, cell_size=8);

  /// [Meshing]
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
  /// [Meshing]

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);
  std::ofstream bin_file("out.binary.cgal", std::ios_base::binary);
  CGAL::Mesh_3::save_binary_file(bin_file, c3t3);

  return 0;
}
