
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Mesh_3/initialize_triangulation_from_gray_image.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <functional>

typedef float Image_word_type;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Parallel tag
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

namespace params = CGAL::parameters;

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("images/skull_2.9.inr");
  /// [Load image]
  CGAL::Image_3 image;
  if (!image.read(fname)) {
    std::cerr << "Error: Cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }
  /// [Domain creation]
  Mesh_domain domain =
    Mesh_domain::create_gray_image_mesh_domain(image, params::iso_value(2.9f).value_outside(0.f));
  /// [Domain creation]

  /// [Mesh criteria]
  Mesh_criteria criteria(params::facet_angle(30).facet_size(6).facet_distance(2).
                                 cell_radius_edge_ratio(3).cell_size(8));

  /// [Meshing]
  C3t3 c3t3;
  initialize_triangulation_from_gray_image(c3t3,
                                           domain,
                                           image,
                                           criteria,
                                           2.9f,//isolevel
                                           Image_word_type(0));
  CGAL::refine_mesh_3(c3t3, domain, criteria);
  /// [Meshing]

  /// Output
  CGAL::dump_c3t3(c3t3, "out");

  return 0;
}
