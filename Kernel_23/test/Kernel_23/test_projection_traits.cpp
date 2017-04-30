#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>

#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Projection_traits_xy_3<Epic>                  K;

typedef K::Line_2                                           Line_2;
typedef K::Point_2                                          Point_2;
typedef K::Weighted_point_2                                 Weighted_point_2;
typedef K::Vector_2                                         Vector_2;

int main()
{
  Point_2 p3(1,0,1), q3(1,1,2), r3(1,2,3), s3(2,1,4);
  Weighted_point_2 wp3(p3, 1), wq3(q3, 2), wr3(r3, 3), ws3(s3, -4);
  Vector_2 v3(1, 0, 0), w3(0,1,0);

  Epic::Point_2 p2(1,0), q2(1,1), r2(1,2), s2(2,1);
  Epic::Weighted_point_2 wp2(p2, 1), wq2(q2, 2), wr2(r2, 3), ws2(s2, -4);
  Epic::Vector_2 v2(1, 0), w2(0,1);

  K k;

  assert( k.compute_scalar_product_2_object()(v3, w3) ==
          v2 * w2 );

  assert( k.collinear_2_object()(p3,q3,r3) ==
          CGAL::collinear(p2,q2,r2) );

  assert( k.collinear_are_ordered_along_line_2_object()(p3,q3,r3) ==
          CGAL::collinear_are_ordered_along_line(p2,q2,r2) );

  assert( k.compute_squared_length_2_object()(v3) == v2.squared_length() );

  assert( p3 == k.construct_point_2_object()(p3) );
  assert( p3 == k.construct_point_2_object()(wp3) );
  assert( wp3 == k.construct_weighted_point_2_object()(p3, 1) );
  assert( wp3 == k.construct_weighted_point_2_object()(wp3) );

  assert( k.compare_power_distance_2_object()(p3, wq3, wr3) ==
           CGAL::compare_power_distance(p2, wq2, wr2) );

  assert( k.compute_power_product_2_object()(wp3, wq3) ==
            CGAL::power_product(wp2, wq2) );

  assert( k.compute_squared_radius_smallest_orthogonal_circle_2_object()(wp3, wq3, wr3) ==
            CGAL::squared_radius_smallest_orthogonal_circle(wp2, wq2, wr2) );
  assert( k.compute_squared_radius_smallest_orthogonal_circle_2_object()(wp3, wq3) ==
            CGAL::squared_radius_smallest_orthogonal_circle(wp2, wq2) );
  assert( k.compute_squared_radius_smallest_orthogonal_circle_2_object()(wp3, wq3) ==
            CGAL::squared_radius_smallest_orthogonal_circle(wp2, wq2) );

  Line_2 l = k.construct_radical_axis_2_object()(wp3, wq3);
  Point_2 p = l.point(42);
  Epic::Point_2 ep(p.x(), p.y());
  assert( CGAL::compare_power_distance(ep, wp2, wq2) == CGAL::EQUAL );

  Point_2 c = k.construct_weighted_circumcenter_2_object()(wp3, wq3, ws3);
  Epic::Point_2 ec = CGAL::weighted_circumcenter(wp2, wq2, ws2);
  assert( c.x() == ec.x() && c.y() == ec.y() );

  assert( k.power_side_of_bounded_power_circle_2_object()(wp3, wq3, ws3, wr3) ==
            CGAL::power_side_of_bounded_power_circle(wp2, wq2, ws2, wr2) );
  assert( k.power_side_of_bounded_power_circle_2_object()(wp3, wq3, ws3) ==
            CGAL::power_side_of_bounded_power_circle(wp2, wq2, ws2) );
  assert( k.power_side_of_bounded_power_circle_2_object()(wp3, wq3) ==
            CGAL::power_side_of_bounded_power_circle(wp2, wq2) );

  assert( k.power_side_of_oriented_power_circle_2_object()(wp3, wq3, ws3, wr3) ==
            CGAL::power_side_of_oriented_power_circle(wp2, wq2, ws2, wr2) );
  assert( k.power_side_of_oriented_power_circle_2_object()(wp3, wq3, wr3) ==
            CGAL::power_side_of_oriented_power_circle(wp2, wq2, wr2) );
  assert( k.power_side_of_oriented_power_circle_2_object()(wp3, wq3) ==
            CGAL::power_side_of_oriented_power_circle(wp2, wq2) );

  std::cout << "done" << std::endl;
 return 0;
}
