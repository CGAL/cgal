#include <vector>
#include <list>

#include <CGAL/config.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>

///////////////
//circle segment traits
//////////////
typedef CGAL::Quotient<CGAL::MP_Float>               Number_type;
typedef CGAL::Cartesian<Number_type>                 Kernel;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>    Sub_traits_2;
typedef CGAL::Arr_polyline_traits_2<Sub_traits_2>    Traits_2;
typedef Sub_traits_2::CoordNT                        CoordNT;
typedef Sub_traits_2::Point_2                        Point_2;
typedef Sub_traits_2::Curve_2                        Subcurve_2;
typedef Sub_traits_2::X_monotone_curve_2             Subcurve_x_monotone_2;
typedef boost::variant<Point_2, Subcurve_x_monotone_2>
  Make_sub_x_monotone_result;

void check_equal()
{
  Traits_2 Polycurve_traits_2;
  auto equal_2 = Polycurve_traits_2.equal_2_object();
  auto construct_x_monotone_curve_2 =
    Polycurve_traits_2.construct_x_monotone_curve_2_object();

  Subcurve_2 curve1, curve2, curve3;

  Kernel::Point_2 p1 = Kernel::Point_2(-5, 0);
  Kernel::Point_2 mid = Kernel::Point_2(0, 5);
  Kernel::Point_2 p2 = Kernel::Point_2(5, 0);
  curve1= Subcurve_2(p1, mid, p2);
  curve2= Subcurve_2(p1, mid, p2);

  // //make x_monotone
  // Traits_2::X_monotone_curve_2 xmc1 =
  //   construct_x_monotone_curve_2(curve1);
  // Traits_2::X_monotone_curve_2 xmc2 =
  //   construct_x_monotone_curve_2(curve2);
  // Traits_2::X_monotone_curve_2 xmc3 =
  //   construct_x_monotone_curve_2(curve3);

  // //simple equal
  // bool Are_equal = equal_2(xmc1, xmc1);
  // std::cout << "Equal_2::Two equal semi circles are computed as: "
  //           << ((Are_equal) ? "Equal" : "Not equal") << std::endl;

  // Are_equal = equal_2(xmc1, xmc3);
  // std::cout << "Equal_2::Two un-equal semi circles are computed as: "
  //           << ((Are_equal) ? "Equal" : "Not equal") << std::endl;
}

void check_intersect(Traits_2::Make_x_monotone_2
                     make_x_monotone_2,
                     Traits_2::Intersect_2 intersect_2)
{
  Subcurve_2 curve1, curve2, curve3;

  // Create a circular arc that correspond to the upper half of the
  // circle centered at (1,1) with squared radius 3. We create the
  // circle with clockwise orientation, so the arc is directed from
  // (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
  Kernel::Point_2 c1 = Kernel::Point_2(1, 1);
  Kernel::Circle_2 circ1 = Kernel::Circle_2(c1, 3, CGAL::CLOCKWISE);
  CoordNT one_minus_sqrt_3 = CoordNT(1, -1, 3);
  CoordNT one_plus_sqrt_3 = CoordNT(1, 1, 3);
  Point_2 s1 = Point_2(one_minus_sqrt_3, CoordNT(1));
  Point_2 t1 = Point_2(one_plus_sqrt_3, CoordNT(1));
  curve1 = Subcurve_2(circ1, s1, t1);
  curve2 = Subcurve_2(circ1, s1, t1);
  //push the same semi circle again
  //curves.push_back(Subcurve_2(circ1, s1, t1));

  //make x_monotone
  std::vector<Make_sub_x_monotone_result> x_monotone_curves;
  make_x_monotone_2(curve1, std::back_inserter(x_monotone_curves));
  make_x_monotone_2(curve2, std::back_inserter(x_monotone_curves));

  const auto* x_curve1 = boost::get<Subcurve_x_monotone_2>(curves[0]);
  const auto* x_curve2 = boost::get<Subcurve_x_monotone_2>(curves[1]);

  std::vector<Intersection_result> Points_of_intersection;

  //intersect_2(X_monotone_curve1, X_monotone_curve2,
  //            std::back_inserter(Points_of_intersection));

 // Create a circular arc of the unit circle, directed clockwise from
  // (-1/2, sqrt(3)/2) to (1/2, sqrt(3)/2). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2 c6 = Kernel::Point_2(0, 0);
  CoordNT sqrt_3_div_2 =
    CoordNT(Number_type(0), Number_type(1,2), Number_type(3));
  Point_2 s6 = Point_2(Number_type(-1, 2), sqrt_3_div_2);
  Point_2 t6 = Point_2(Number_type(1, 2), sqrt_3_div_2);

  curve3 = Subcurve_2(c6, 1, CGAL::CLOCKWISE, s6, t6);
  make_x_monotone_2(curve3, std::back_inserter(x_curves));
  const auto* x_curve2 = boost::get<>(&x_curves[2]);

  Points_of_intersection.clear();
  //intersect_2(X_monotone_curve1, X_monotone_curve2,
  //            std::back_inserter(Points_of_intersection));
}

void
check_compare_end_points_xy_2(Traits_2::Compare_endpoints_xy_2
                              compare_endpoints_xy_2,
                              Traits_2::Make_x_monotone_2
                              make_x_monotone_2)
{
  Subcurve_2 curve1, curve2;

  // Create a circular arc that correspond to the upper half of the
  // circle centered at (1,1) with squared radius 3. We create the
  // circle with clockwise orientation, so the arc is directed from
  // (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
  Kernel::Point_2 c1 = Kernel::Point_2(1, 1);
  Kernel::Circle_2 circ1 = Kernel::Circle_2(c1, 3, CGAL::CLOCKWISE);
  CoordNT one_minus_sqrt_3 = CoordNT(1, -1, 3);
  CoordNT one_plus_sqrt_3 = CoordNT(1, 1, 3);
  Point_2 s1 = Point_2(one_minus_sqrt_3, CoordNT(1));
  Point_2 t1 = Point_2(one_plus_sqrt_3, CoordNT(1));
  curve1= Subcurve_2(circ1, s1, t1);

  //make x_monotone
  std::vector<Make_sub_x_monotone_result> x_curves;
  make_x_monotone_2(curve1, std::back_inserter(x_curves));

  const auto* x_curve1 = boost::get<Subcurve_x_monotone_2>(&x_curves[0]);
  const auto* x_curve2 = boost::get<Subcurve_x_monotone_2>(&x_curves[1]);

  auto res = compare_endpoints_xy_2(x_curve1);
  std::cout<< "The first result is: " << res << std::endl;

  Kernel::Point_2 c2 = Kernel::Point_2(1, 1);
  Kernel::Circle_2 circ2 = Kernel::Circle_2(c2, 3, CGAL::COUNTERCLOCKWISE);
  Point_2 t2 = Point_2(one_minus_sqrt_3, CoordNT(1));
  Point_2 s2 = Point_2(one_plus_sqrt_3, CoordNT(1));
  curve2 = Subcurve_2(circ2, s1, t1);
  make_x_monotone_2(curve2, std::back_inserter(x_curves));
  const auto* x_curve2 = boost::get<Subcurve_x_monotone_2>(&x_curves[1]);

  res = compare_endpoints_xy_2(x_curve2);

  std::cout<< "The second result is: " << res << std::endl;
}

void check_split(Traits_2::Split_2  split_2,
                 Traits_2::Make_x_monotone_2 make_x_monotone_2)
{
  Subcurve_2 curve;

  // Create a circular arc that correspond to the upper half of the
  // circle centered at (1,1) with squared radius 3. We create the
  // circle with clockwise orientation, so the arc is directed from
  // (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
  Kernel::Point_2 c1 = Kernel::Point_2(1, 1);
  Kernel::Circle_2 circ1 = Kernel::Circle_2(c1, 3, CGAL::CLOCKWISE);
  CoordNT one_minus_sqrt_3 = CoordNT(1, -1, 3);
  CoordNT one_plus_sqrt_3 = CoordNT(1, 1, 3);
  Point_2 s1 = Point_2(one_minus_sqrt_3, CoordNT(1));
  Point_2 t1 = Point_2(one_plus_sqrt_3, CoordNT(1));
  curve= Subcurve_2(circ1, s1, t1);

  //make x_monotone
  std::vector<Make_sub_x_monotone_result> x_curves;
  make_x_monotone_2(curve, std::back_inserter(x_curves));
  const auto* x_curve1 = boost::get<Subcurve_x_monotone_2>(&x_curves[0]);

  // Subcurve_x_monotone_2 split_x_monotone_curve1, split_x_monotone_curve2 ;
  //split_2(X_monotone_curve, Kernel::Point_2::Kernel::Point_2(1, 4),
  //        split_x_monotone_curve1, split_x_monotone_curve2);
}

void check_is_vertical(Traits_2::Make_x_monotone_2
                       make_x_monotone_2,
                       Traits_2::Is_vertical_2 is_vertical)
{
  std::vector<Subcurve_2> curves;

  // Create a circular arc defined by two endpoints and a midpoint,
  // all having rational coordinates. This arc is the upper-right
  // quarter of a circle centered at the origin with radius 5.
  Kernel::Point_2 p1 = Kernel::Point_2(0, 5);
  Kernel::Point_2 mid = Kernel::Point_2(3, 4);
  Kernel::Point_2 p2 = Kernel::Point_2(5, 0);
  Kernel::Point_2 p3 = Kernel::Point_2(0, -5);
  curves.push_back(Subcurve_2(p1, mid, p2)); //quarter of a circle
  curves.push_back(Subcurve_2(p1, mid, p3));  //semi-circle

  //convert all curves to x-monotone curves
  std::vector<Make_sub_x_monotone_result> x_curves;
  for (const auto& cv : curves)
    make_x_monotone_2(cv, std::back_inserter(x_curves));

  //std::vector<Subcurve_x_monotone_2> x_monotone_polycurves;

  const auto* x_polycurve1 = boost::get<Subcurve_x_monotone_2>(&x_curves[0]);
  const auto* x_polycurve2 = boost::get<Subcurve_x_monotone_2>(&x_curves[1]);

  bool res = is_vertical(x_polycurve1);
  std::cout << "Is_verticle:: The xmonotone curve (quarter circle) is : "
            << ((res)? "vertical" : "not vertical") << std::endl;

  res = is_vertical(x_polycurve2);
  std::cout << "Is_verticle:: The xmonotone curve (Smi-circle) is : "
            << ((res)? "vertical" : "not vertical") << std::endl;
}

void check_compare_y_at_x_2(Traits_2::Make_x_monotone_2 make_x_monotone_2,
                            Traits_2::Compare_y_at_x_2 cmp_y_at_x_2)
{
  std::vector<Subcurve_2> curves;

  // Create a circular arc defined by two endpoints and a midpoint,
  // all having rational coordinates. This arc is the upper-right
  // quarter of a circle centered at the origin with radius 5.
  Kernel::Point_2 p1 = Kernel::Point_2(1, 1);
  Kernel::Point_2 mid = Kernel::Point_2(4, 4);
  Kernel::Point_2 p2 = Kernel::Point_2(7, 1);
  Kernel::Point_2 p3 = Kernel::Point_2(1, 4);
  curves.push_back(Subcurve_2(p1, mid, p2)); //quarter of a circle
  curves.push_back(Subcurve_2(p1, mid, p3));  //semi-circle

  //convert all curves to x-monotone curves
  std::vector<Make_sub_x_monotone_result> x_curves;
  for (const auto& cv : curves)
    make_x_monotone_2(cv, std::back_inserter(x_curves));

  const auto* x_polycurve1 = boost::get<Subcurve_x_monotone_2>(&x_curves[0]);
  const auto* x_polycurve2 = boost::get<Subcurve_x_monotone_2>(&x_curves[1]);

  Kernel::Point_2 p_test = Kernel::Point_2(3, 1);

  // int res = cmp_y_at_x_2(p_test, x_monotone_polycurve1);
  // cmp_y_at_x_2(x_monotone_polycurve1, CGAL::ARR_MIN_END,
  //              x_monotone_polycurve2);
}

void check_push_back(Traits_2::Make_x_monotone_2
                     make_x_monotone_2,
                     Traits_2::Push_back_2  push_back_2)
{
  std::vector<Subcurve_2> curves;

  //check if segment is pushed in empty curve.
  Kernel::Point_2 p1 = Kernel::Point_2(1, 1);
  Kernel::Point_2 mid = Kernel::Point_2(4, 4);
  Kernel::Point_2 p2 = Kernel::Point_2(7, 1);

  Kernel::Point_2 mid2 = Kernel::Point_2(10, 3);
  Kernel::Point_2 p3 = Kernel::Point_2(7, 7);

  curves.push_back(Subcurve_2(p1, mid, p2));
  curves.push_back(Subcurve_2(p2, mid2, p3));

  CGAL::internal::Polycurve_2<Subcurve_2, Point_2> polycurve;

  ////pushing segments in polycurve
  push_back_2(polycurve, curves[0]);
  std::cout << "size of polycurve after 1 push_back: "
            << polycurve.number_of_subcurves() << std::endl;

  push_back_2(polycurve, curves[1]);
  //throws a warning "size is depricated"
  std::cout << "size of polycurve after 2 push_backs: "
            << polycurve.number_of_subcurves() << std::endl;
}

int main()
{
  Traits_2 Polycurve_traits_2;

  // Compare_x_2
  //  Traits_2::Compare_x_2 compare_x_2 =
  //   Polycurve_traits_2.compare_xy_2_object();

  // number of points
  // Traits_2::Number_of_points_2 num_of_points =
  //   Polycurve_traits_2.number_of_points_2_object();

  //construct min vertex
  auto cnst_min_vertex = Polycurve_traits_2.construct_min_vertex_2_object();

  //construct max vertex
  auto cnst_max_vertex_2 = Polycurve_traits_2.construct_max_vertex_2_object();

  //is vertical (return bool)
  auto is_vertical = Polycurve_traits_2.is_vertical_2_object();

  //Compare y at x 2 (return comparison_result)
  auto cmp_y_at_x_2 = Polycurve_traits_2.compare_y_at_x_2_object();

  //compare y at x left
  auto cmp_y_at_x_left_2 = Polycurve_traits_2.compare_y_at_x_left_2_object();

  //compare y at x right
  auto cmp_y_at_x_right_2 = Polycurve_traits_2.compare_y_at_x_right_2_object();

  //equal_2
  auto equal_2 = Polycurve_traits_2.equal_2_object();

  //compare end points xy_2
  auto compare_endpoints_xy_2 =
    Polycurve_traits_2.compare_endpoints_xy_2_object();

  //construct opposite
  auto construct_opposite_2 = Polycurve_traits_2.construct_opposite_2_object();

  //make x_monotone
  auto make_x_monotone_2 = Polycurve_traits_2.make_x_monotone_2_object();

  //push back
  auto push_back_2 = Polycurve_traits_2.push_back_2_object();

  //push front
  auto push_front_2 = Polycurve_traits_2.push_front_2_object();

  //split_2
  auto split_2 = Polycurve_traits_2.split_2_object();

  //Intersect_2
  auto intersect_2 = Polycurve_traits_2.intersect_2_object();

  //Are_mergable
  auto are_mergeable_2 = Polycurve_traits_2.are_mergeable_2_object();

  //Merge_2
  auto merge_2 = Polycurve_traits_2.merge_2_object();

  //construct_curve_2
  auto construct_curve_2 = Polycurve_traits_2.construct_curve_2_object();

  //construct x_monotone curve_2
  auto construct_x_monotone_curve_2 =
    Polycurve_traits_2.construct_x_monotone_curve_2_object();

  //check_equal();
  check_intersect(make_x_monotone_2, intersect_2);
  check_compare_end_points_xy_2(compare_endpoints_xy_2, make_x_monotone_2);
  check_split(split_2, make_x_monotone_2);
  check_is_vertical(make_x_monotone_2, is_vertical);
  check_compare_y_at_x_2(make_x_monotone_2, cmp_y_at_x_2);
  check_push_back(make_x_monotone_2, push_back_2);

  return 0;
}
