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

// Custom_initial_points_generator will put points on the mesh for initialisation.
// Those points are objects of type std::tuple<Weighted_point_3, int, Index>.
//   Weighted_point_3 is the point's position and weight,
//   int is the dimension of the minimal dimension subcomplex on which the point lies,
//   Index is the underlying subcomplex index.
struct Custom_initial_points_generator
{
  const CGAL::Image_3& image_;

  template <typename OutputIterator>
  OutputIterator operator()(OutputIterator pts) const
  {
    typedef Tr::Geom_traits     Gt;
    typedef Gt::Point_3         Point_3;
    typedef Gt::Vector_3        Vector_3;
    typedef Gt::Segment_3       Segment_3;
    typedef Mesh_domain::Index  Index;

    typename C3t3::Triangulation::Geom_traits::Construct_weighted_point_3 cwp =
        c3t3.triangulation().geom_traits().construct_weighted_point_3_object();

    // Add points along the segment
    Segment_3 segment(Point_3(  0.0, 50.0, 66.66),
                      Point_3(100.0, 50.0, 66.66));


    Point_3 source = segment.source();
    Vector_3 vector = segment.to_vector();
    double edge_size = 5;
    std::size_t nb = static_cast<int>(CGAL::sqrt(segment.squared_length()) / edge_size);
    const double frac = 1. / (double)nb;

    for (std::size_t i = 1; i < nb; i++)
    {
      *pts++ = {cwp(source + (i * frac) * vector), 1, Index(1));
    }
    return pts;
  }
};

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
      , params::features_detector(CGAL::Mesh_3::Detect_features_on_image_bbox())
  );

  // Mesh criteria
  Mesh_criteria criteria(params::facet_angle(30).facet_size(3).facet_distance(1).edge_size(3)
                         .cell_radius_edge_ratio(3).cell_size(3)
  );

  /// [Meshing]
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria
    , params::initial_points_generator(Custom_initial_points_generator{ image })
  );
  /// [Meshing]

  // Output
  CGAL::dump_c3t3(c3t3, "out");

  return 0;
}
