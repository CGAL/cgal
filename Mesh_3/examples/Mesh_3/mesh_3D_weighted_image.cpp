#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/generate_label_weights.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/IO/File_binary_mesh_3.h>
#include <CGAL/tags.h>

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh_domain = CGAL::Labeled_mesh_domain_3<K>;

// Triangulation
using Tr   = CGAL::Mesh_triangulation_3<Mesh_domain,
                                        CGAL::Default,
                                        CGAL::Parallel_if_available_tag>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char* argv[])
{
  /// [Loads image]
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("images/liver.inr.gz");
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  /// [Loads image]

  /// [Domain creation]
  const float sigma = 10.f;
  CGAL::Image_3 img_weights =
    CGAL::Mesh_3::generate_label_weights(image, sigma);

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
  CGAL::IO::save_binary_file(bin_file, c3t3);

  return 0;
}
