#include "random_labeled_image.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Mesh_3/Construct_initial_points_labeled_image.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

#include <CGAL/SMDS_3/Dump_c3t3.h>

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

namespace params = CGAL::parameters;

template<class C3T3, class MeshDomain>
void initialize_triangulation_from_labeled_image(C3T3& c3t3,
      const MeshDomain&   domain,
      const CGAL::Image_3& image)
{
  typedef typename C3T3::Triangulation       Tr;
  typedef typename Tr::Geom_traits           GT;
  typedef typename Tr::Weighted_point        Weighted_point;
  typedef typename Tr::Vertex_handle         Vertex_handle;
  typedef typename MeshDomain::Point_3       Point_3;
  typedef typename MeshDomain::Index         Index;

  typedef typename std::pair<Point_3, Index> ConstructedPoint;

  Tr& tr = c3t3.triangulation();

  typename GT::Construct_weighted_point_3 cwp =
    tr.geom_traits().construct_weighted_point_3_object();

  std::vector<ConstructedPoint> constructedPoints;

  CGAL::Construct_initial_points_labeled_image construct(image);
  construct(std::back_inserter(constructedPoints), domain, c3t3);

  std::cout << "  " << constructedPoints.size() << " constructed points" << std::endl;

  for (const ConstructedPoint & constructedPoint : constructedPoints)
  {
    const Point_3& point = constructedPoint.first;
    const Index&   index = constructedPoint.second;

    Weighted_point pi = cwp(point);

    /// The following lines show how to insert initial points in the
    /// `c3t3` object. [insert initial points]
    Vertex_handle v = tr.insert(pi);
    // `v` could be null if `pi` is hidden by other vertices of `tr`.
    CGAL_assertion(v != Vertex_handle());
    c3t3.set_dimension(v, 2); // by construction, points are on surface
    c3t3.set_index(v, index);
    /// [insert initial points]
  }
}

int main()
{
  /// [Create the image]
  CGAL::Image_3 image = random_labeled_image();
  /// [Create the image]

  // Domain
  Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image);

  // Mesh criteria
  Mesh_criteria criteria(params::facet_angle(30).facet_size(3).facet_distance(1).
                         cell_radius_edge_ratio(3).cell_size(3));

  /// [Meshing]
  C3t3 c3t3;
  initialize_triangulation_from_labeled_image(c3t3,
                                              domain,
                                              image);
  CGAL::refine_mesh_3(c3t3, domain, criteria);
  /// [Meshing]

  // Output
  CGAL::dump_c3t3(c3t3, "out");

  return 0;
}
