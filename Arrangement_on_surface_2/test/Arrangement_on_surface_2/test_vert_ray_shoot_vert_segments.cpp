#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_simple_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>

template <typename VerticalRayShoot>
void vertical_ray_shooting_query
(const VerticalRayShoot& vrs, bool shoot_up,
 const typename VerticalRayShoot::Arrangement_2::Point_2& q)
{
  // Perform the point-location query.
  auto obj = (shoot_up) ? vrs.ray_shoot_up(q) : vrs.ray_shoot_down(q);

  // Print the result.
  typedef typename VerticalRayShoot::Arrangement_2      Arrangement_2;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

  if (shoot_up) std::cout << "  Shooting up from (" << q << ") : ";
  else std::cout << "  Shooting down from (" << q << ") : ";

  if (auto* e_p = std::get_if<Halfedge_const_handle>(&obj)) {
    // We hit an edge:
    std::cout << "hit an edge: " << (*e_p)->curve() << std::endl;
  }
  else if (auto* v_p = std::get_if<Vertex_const_handle>(&obj)) {
    // We hit a vertex:
    if ((*v_p)->is_isolated())
      std::cout << "hit an isolated vertex: " << (*v_p)->point() << std::endl;
    else std::cout << "hit a vertex: " << (*v_p)->point() << std::endl;
  }
  else if (auto* f_p = std::get_if<Face_const_handle>(&obj)) {
    // We did not hit anything:
    assert((*f_p)->is_unbounded());

    std::cout << "hit nothing." << std::endl;
  }
  else CGAL_error_msg( "Invalid object.");
}

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::FT                                              Number_type;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef Traits_2::X_monotone_curve_2                            Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;
typedef CGAL::Arr_simple_point_location<Arrangement_2>          Simple_pl;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;
typedef CGAL::Arr_trapezoid_ric_point_location<Arrangement_2>   RIC_pl;

int main()
{
  // Construct the arrangement.
  Arrangement_2 arr;
  Simple_pl simple_pl(arr);
  Walk_pl walk_pl(arr);
  RIC_pl ric_pl(arr);

  insert(arr, Segment_2(Point_2(0, 0), Point_2(0, 1)));
  insert(arr, Segment_2(Point_2(0, 3), Point_2(0, 4)));
  insert(arr, Segment_2(Point_2(-1, -1), Point_2(1, -1)));
  insert(arr, Segment_2(Point_2(-1, 5), Point_2(1, 5)));

  // Perform the vertical ray-shooting queries.
  std::cout << "First round of queries: " << std::endl;

  Point_2 q1(0, 2);

  vertical_ray_shooting_query(simple_pl, true, q1);
  vertical_ray_shooting_query(walk_pl, true, q1);
  vertical_ray_shooting_query(ric_pl, true, q1);

  vertical_ray_shooting_query(simple_pl, false, q1);
  vertical_ray_shooting_query(walk_pl, false, q1);
  vertical_ray_shooting_query(ric_pl, false, q1);

  Point_2 q2(0, 1);

  vertical_ray_shooting_query(simple_pl, true, q2);
  vertical_ray_shooting_query(walk_pl, true, q2);
  vertical_ray_shooting_query(ric_pl, true, q2);

  vertical_ray_shooting_query(simple_pl, false, q2);
  vertical_ray_shooting_query(walk_pl, false, q2);
  vertical_ray_shooting_query(ric_pl, false, q2);

  Point_2 q3(0, 3);

  vertical_ray_shooting_query(simple_pl, true, q3);
  vertical_ray_shooting_query(walk_pl, true, q3);
  vertical_ray_shooting_query(ric_pl, true, q3);

  vertical_ray_shooting_query(simple_pl, false, q3);
  vertical_ray_shooting_query(walk_pl, false, q3);
  vertical_ray_shooting_query(ric_pl, false, q3);

  // Insert additional curves and perform the ray-shooting queries again.
  insert (arr, Segment_2(Point_2(-1, 0), Point_2(1, 0)));
  insert (arr, Segment_2(Point_2(-1, 4), Point_2(1, 4)));

  std::cout << "Second round of queries: " << std::endl;

  vertical_ray_shooting_query(simple_pl, true, q1);
  vertical_ray_shooting_query(walk_pl, true, q1);
  vertical_ray_shooting_query(ric_pl, true, q1);

  vertical_ray_shooting_query(simple_pl, false, q1);
  vertical_ray_shooting_query(walk_pl, false, q1);
  vertical_ray_shooting_query(ric_pl, false, q1);

  vertical_ray_shooting_query(simple_pl, true, q2);
  vertical_ray_shooting_query(walk_pl, true, q2);
  vertical_ray_shooting_query(ric_pl, true, q2);

  vertical_ray_shooting_query(simple_pl, false, q2);
  vertical_ray_shooting_query(walk_pl, false, q2);
  vertical_ray_shooting_query(ric_pl, false, q2);

  vertical_ray_shooting_query(simple_pl, true, q3);
  vertical_ray_shooting_query(walk_pl, true, q3);
  vertical_ray_shooting_query(ric_pl, true, q3);

  vertical_ray_shooting_query(simple_pl, false, q3);
  vertical_ray_shooting_query(walk_pl, false, q3);
  vertical_ray_shooting_query(ric_pl, false, q3);

  return (0);
}
