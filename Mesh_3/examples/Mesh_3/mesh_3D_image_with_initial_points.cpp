#include "random_labeled_image.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

#include <CGAL/Mesh_3/Detect_features_in_image.h>

#include <CGAL/IO/File_medit.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Image_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Image_domain> Mesh_domain;

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

namespace params = CGAL::parameters;

int main()
{
  const std::string fname = CGAL::data_file_path("images/420.inr");
  // Loads image
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }

  // Domain
  Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image
      , params::features_detector(CGAL::Mesh_3::Detect_features_in_image()));

  // Mesh criteria
  Mesh_criteria criteria(params::facet_angle(30).facet_size(3).facet_distance(1).edge_size(3)
                         .cell_radius_edge_ratio(3).cell_size(3));

  using Point_3 = K::Point_3;
  using Weighted_point_3 = K::Weighted_point_3;
  using Index = Mesh_domain::Index;
  using Initial_point_t = std::tuple<Weighted_point_3, int, Index>;

  // Creation of the initial_points vector
  std::vector<Initial_point_t> initial_points = {
    {Weighted_point_3(Point_3(30.0, 50.0, 83.33), 30.0), 1, Index(1)},
    {Weighted_point_3(Point_3(70.0, 50.0, 83.33), 50.0), 1, Index(1)}
  };

  /// [Meshing]
  // Mesh generation from labeled image with initial points.
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      params::initial_points(std::cref(initial_points)));
  /// [Meshing]

  // Output
  std::ofstream ofs("out.mesh");
  CGAL::IO::write_MEDIT(ofs, c3t3);

  return 0;
}
