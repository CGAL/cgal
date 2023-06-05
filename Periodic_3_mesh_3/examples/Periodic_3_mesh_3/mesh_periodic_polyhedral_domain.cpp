#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {

// This is a wrapper to convert a polyhedral surface to a periodic polyhedral domain
// over a user-provided canonical domain.
//
// It is the user's responsability to ensure that the polyhedral domain is actually periodic
// over the canonical domain, i.e. there is periodic continuity at the boundaries
// of the canonical domain.
template<class TriangleMesh, class K>
class Polyhedral_to_periodic_labeling_function_wrapper
{
public:
  using return_type = int;
  using Point_3 = typename K::Point_3;
  using GT = typename details::Periodic_3_mesh_geom_traits_generator<K>::type;

private:
  const TriangleMesh& m_tmesh;
  CGAL::Side_of_triangle_mesh<TriangleMesh, GT> m_sotm;
  GT m_gt;

public:
  explicit Polyhedral_to_periodic_labeling_function_wrapper(const TriangleMesh& tmesh,
                                                            const CGAL::Iso_cuboid_3<K>& domain)
    : m_tmesh(tmesh), m_sotm(m_tmesh), m_gt(domain)
  {
    CGAL_precondition(CGAL::is_closed(tmesh));
  }

  Polyhedral_to_periodic_labeling_function_wrapper(const Polyhedral_to_periodic_labeling_function_wrapper& other)
    : m_tmesh(other.m_tmesh), m_sotm(m_tmesh), m_gt(other.m_gt)
  { }

  return_type operator()(const Point_3& p) const
  {
    const Point_3 cp = P3T3::internal::robust_canonicalize_point(p, m_gt);
    CGAL::Bounded_side res = m_sotm(cp);
    return ((res == ON_BOUNDED_SIDE) ? 1 : 2); // set a region to '0' if it is not to be meshed
  }
};

} // namespace CGAL

// Kernel
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point = K::Point_3;
using Iso_cuboid = K::Iso_cuboid_3;

// Domain
using Mesh = CGAL::Surface_mesh<Point>;
using Periodic_polyhedral_domain = CGAL::Polyhedral_to_periodic_labeling_function_wrapper<Mesh, K>;

// Optional: add polyline features to be protected (i.e., preserved in the output)
// using Periodic_mesh_domain = CGAL::Labeled_mesh_domain_3<K>; // no feature protection
using Periodic_mesh_domain = CGAL::Labeled_mesh_domain_3<K>;
using Polyline_3 = std::vector<Point>;
using Polylines = std::list<Polyline_3>;

// Triangulation
using Tr = CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Periodic_mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

// To avoid verbose function and named parameters call
namespace params = CGAL::parameters;

// An arbitrary, simple polyhedral shape
void generate_periodic_diamond(const Point& origin, // bottom, front, left point of the canonical domain
                               const FT xs, const FT ys, const FT zs,
                               Mesh& sm)
{
  using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;

  auto vpm = get(CGAL::vertex_point, sm);

  // Three points on in the middle of the canonical domain
  vertex_descriptor v0 = add_vertex(sm);
  vertex_descriptor v1 = add_vertex(sm);
  vertex_descriptor v2 = add_vertex(sm);
  put(vpm, v0, Point(origin.x() + 0.25 * xs, origin.y() + 0.25 * ys, origin.z() + 0.5 * zs));
  put(vpm, v1, Point(origin.x() + 0.75 * xs, origin.y() + 0.25 * ys, origin.z() + 0.5 * zs));
  put(vpm, v2, Point(origin.x() + 0.50 * xs, origin.y() + 0.75 * ys, origin.z() + 0.5 * zs));

  // two points out of the domain, but the intersection of this diamond
  // with the canonical domain forms a periodic domain
  vertex_descriptor v_up = add_vertex(sm);
  vertex_descriptor v_do = add_vertex(sm);
  put(vpm, v_up, Point(origin.x() + 0.5 * xs, origin.y() + 0.5 * ys, origin.z() - 0.5 * zs));
  put(vpm, v_do, Point(origin.x() + 0.5 * xs, origin.y() + 0.5 * ys, origin.z() + 1.5 * zs));

  CGAL::Euler::add_face(std::initializer_list<vertex_descriptor>{v1, v_up, v0}, sm);
  CGAL::Euler::add_face(std::initializer_list<vertex_descriptor>{v2, v_up, v1}, sm);
  CGAL::Euler::add_face(std::initializer_list<vertex_descriptor>{v0, v_up, v2}, sm);

  CGAL::Euler::add_face(std::initializer_list<vertex_descriptor>{v0, v_do, v1}, sm);
  CGAL::Euler::add_face(std::initializer_list<vertex_descriptor>{v1, v_do, v2}, sm);
  CGAL::Euler::add_face(std::initializer_list<vertex_descriptor>{v2, v_do, v0}, sm);

  CGAL::IO::write_OFF("periodic_spike.off", sm, CGAL::parameters::stream_precision(17));
  assert(is_valid_polygon_mesh(sm));
}

int main(int argc, char** argv)
{
  const FT x_span = (argc > 1) ? atof(argv[1]) : 1;
  const FT y_span = (argc > 2) ? atof(argv[2]) : x_span;
  const FT z_span = (argc > 3) ? atof(argv[3]) : x_span;
  const FT min_span = (std::min)({x_span, y_span, z_span});

  const int number_of_copies_in_output = (argc > 4) ? atoi(argv[4]) : 4; // can be 1, 2, 4, or 8

  Mesh sm;
  generate_periodic_diamond(CGAL::ORIGIN, x_span, y_span, z_span, sm);
  Iso_cuboid canonical_cube(0, 0, 0, x_span, y_span, z_span);

  Periodic_polyhedral_domain ppd(sm, canonical_cube);
  Periodic_mesh_domain domain(ppd, canonical_cube);

  Periodic_mesh_criteria criteria(params::edge_size(min_span)
                                         .facet_angle(30)
                                         .facet_size(0.035 * min_span)
                                         .facet_distance(0.025 * min_span)
                                         .cell_radius_edge_ratio(2.)
                                         .cell_size(0.05 * min_span));

  // Mesh generation
  C3t3 c3t3 = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria);

  std::ofstream medit_file("output_periodic_polyhedral_shape.mesh");
  CGAL::IO::output_periodic_mesh_to_medit(medit_file, c3t3, number_of_copies_in_output);

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}
