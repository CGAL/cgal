#include "random_labeled_image.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

#include <CGAL/SMDS_3/Dump_c3t3.h>


#include <CGAL/Mesh_3/Detect_features_on_image_bbox.h>

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Image_domain = CGAL::Labeled_mesh_domain_3<K>;
using Mesh_domain = CGAL::Mesh_domain_with_polyline_features_3<Image_domain>;

#ifdef CGAL_CONCURRENT_MESH_3
using Concurrency_tag = CGAL::Parallel_tag;
#else
using Concurrency_tag = CGAL::Sequential_tag;
#endif

// Triangulation
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

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
      , params::features_detector(CGAL::Mesh_3::Detect_features_on_image_bbox()));

  // Mesh criteria
  const double edge_size = 3;
  Mesh_criteria criteria(params::edge_size(edge_size)
                         .facet_angle(30).facet_size(3).facet_distance(1)
                         .cell_radius_edge_ratio(3).cell_size(3));

  // custom_initial_points_generator will put points on the mesh for initialization.
  // Those points are objects of type std::tuple<Weighted_point_3, int, Index>.
  //   Weighted_point_3 is the point's position and weight,
  //   int is the dimension of the minimal dimension subcomplex on which the point lies,
  //   Index is the underlying subcomplex index.
  auto custom_initial_points_generator = [&](auto pts_out_iterator, int) {
    using Point_3 = K::Point_3;
    using Weighted_point_3 = K::Weighted_point_3;
    using Index = Mesh_domain::Index;
    using Point_dim_index = std::tuple<Weighted_point_3, int, Index>;

    Point_3 a{0.0, 50.0, 66.66};
    Point_3 b{100.0, 50.0, 66.66};

    // Add points along the segment [a, b]
    double dist_ab = CGAL::sqrt(CGAL::squared_distance(a, b));
    int nb = static_cast<int>(std::floor(dist_ab / edge_size));
    auto vector = (b - a) / static_cast<double>(nb);

    Point_3 p = a;
    for(int i = 0; i < nb; ++i) {
      *pts_out_iterator++ = Point_dim_index{Weighted_point_3{p}, 1, Index(1)};
      p += vector;
    }
    return pts_out_iterator;
  };

  /// [Meshing]
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
      params::initial_points_generator(custom_initial_points_generator));
  /// [Meshing]

  // Output
  CGAL::dump_c3t3(c3t3, "out");

  return 0;
}
