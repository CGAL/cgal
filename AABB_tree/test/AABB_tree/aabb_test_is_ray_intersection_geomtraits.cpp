#include <CGAL/assertions.h>
#include <CGAL/internal/AABB_tree/Is_ray_intersection_geomtraits.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

struct nope {};

struct AABBGeomTraits {
  typedef nope Sphere_3;
  typedef nope Point_3;
  typedef nope Do_intersect_3;
  typedef nope Intersect_3;
  typedef nope Construct_sphere_3;
  typedef nope Compute_closest_point_3;
  typedef nope Has_on_bounded_side_3;
  typedef nope Compute_squared_radius_3;
  typedef nope Compute_squared_distance_3;
  Do_intersect_3 do_intersect_3_object();
  Intersect_3 intersect_3_object();
  Construct_sphere_3 construct_sphere_3_object();
  Compute_closest_point_3 compute_closest_point_3_object();
  Has_on_bounded_side_3 has_on_bounded_side_3_object();
  Compute_squared_radius_3 compute_squared_radius_3_object();
  Compute_squared_distance_3 compute_squared_distance_3_object();
}; /* end AABBGeomTraits */

int main()
{
  using namespace CGAL::internal::AABB_tree;

  CGAL_static_assertion_msg(
    (Is_ray_intersection_geomtraits<CGAL::Epeck>::value),
    "CGAL::Epeck should be a RayIntersectionGeomTraits");
  CGAL_static_assertion_msg(
    (Is_ray_intersection_geomtraits< CGAL::Simple_cartesian<double> >::value),
    "CGAL::Epeck should be a RayIntersectionGeomTraits");
  CGAL_static_assertion_msg(
    (!Is_ray_intersection_geomtraits<AABBGeomTraits>::value),
    "Pure AABBGeomTraits shouldn't be a RayIntersectionGeomTraits");
  CGAL_static_assertion_msg(
    (!Is_ray_intersection_geomtraits<nope>::value),
    "The empty struct shouldn't be a RayIntersectionGeomTraits");

  return 0;
}
