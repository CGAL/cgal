
// Constructing an arrangement of polycurves.

#include <CGAL/basic.h>
#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl; 
  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <vector>
#include <list>

#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>


///////////////
//circle segment traits
//////////////
typedef CGAL::Quotient<CGAL::MP_Float>                                    Number_type;
typedef CGAL::Cartesian<Number_type>                                      Kernel;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>                         Arc_traits_2;
typedef CGAL::Arr_polyline_traits_2<Arc_traits_2>                         Polycurve_arc_traits_2;
typedef Arc_traits_2::CoordNT                                             CoordNT;
typedef Arc_traits_2::Point_2                                             Arc_point_2;
typedef Arc_traits_2::Curve_2                                             Arc_section_2;
typedef Arc_traits_2::X_monotone_curve_2                                  Arc_section_x_monotone_2;
typedef CGAL::Arrangement_2<Polycurve_arc_traits_2>                       Arc_arrangment_2;

void check_equal(Polycurve_arc_traits_2::Equal_2 equal_2, Polycurve_arc_traits_2::Make_x_monotone_2  make_x_monotone_2)
{

  Arc_section_2 curve1, curve2, curve3;

  // Create a circular arc that correspond to the upper half of the
  // circle centered at (1,1) with squared radius 3. We create the
  // circle with clockwise orientation, so the arc is directed from
  // (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
  Kernel::Point_2   c1 = Kernel::Point_2 (1, 1);
  Kernel::Circle_2  circ1 = Kernel::Circle_2 (c1, 3, CGAL::CLOCKWISE);
  CoordNT           one_minus_sqrt_3 = CoordNT (1, -1, 3);
  CoordNT           one_plus_sqrt_3 = CoordNT (1, 1, 3);
  Arc_point_2       s1 = Arc_point_2 (one_minus_sqrt_3, CoordNT (1));
  Arc_point_2       t1 = Arc_point_2 (one_plus_sqrt_3, CoordNT (1));
  curve1= Arc_section_2 (circ1, s1, t1);
  curve2= Arc_section_2 (circ1, s1, t1);
  //push the same semi circle again
  //curves.push_back (Arc_section_2 (circ1, s1, t1));

  //make x_monotone
  std::vector<CGAL::Object> X_monotone_curves;
  make_x_monotone_2(curve1, std::back_inserter(X_monotone_curves));
  make_x_monotone_2(curve2, std::back_inserter(X_monotone_curves));

  Arc_section_x_monotone_2 X_monotone_curve1, X_monotone_curve2, X_monotone_curve3 ;
  CGAL::assign( X_monotone_curve1, X_monotone_curves[ 0 ] );
  CGAL::assign( X_monotone_curve2, X_monotone_curves[ 1 ] );

  //simple equal
  bool Are_equal = equal_2(X_monotone_curve1, X_monotone_curve2);

  std::cout << "The two semi circles are " << ( (Are_equal) ? "equal" : "Not equal") << std::endl;

 // Create a circular arc of the unit circle, directed clockwise from
  // (-1/2, sqrt(3)/2) to (1/2, sqrt(3)/2). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2  c6 = Kernel::Point_2 (0, 0);
  CoordNT          sqrt_3_div_2 = CoordNT (Number_type(0), Number_type(1,2), Number_type(3));
  Arc_point_2          s6 = Arc_point_2 (Number_type (-1, 2), sqrt_3_div_2);
  Arc_point_2          t6 = Arc_point_2 (Number_type (1, 2), sqrt_3_div_2);
  
  curve3 = Arc_section_2 (c6, 1, CGAL::CLOCKWISE, s6, t6);
  make_x_monotone_2(curve3, std::back_inserter(X_monotone_curves));
  CGAL::assign( X_monotone_curve2, X_monotone_curves[ 2 ] );

  Are_equal = equal_2(X_monotone_curve1, X_monotone_curve3);

  std::cout << "The two semi circles are " << ( (Are_equal) ? "equal" : "Not equal") << std::endl;
}

void check_intersect(Polycurve_arc_traits_2::Make_x_monotone_2  make_x_monotone_2,  Polycurve_arc_traits_2::Intersect_2  intersect_2)
{
  Arc_section_2 curve1, curve2, curve3;

  // Create a circular arc that correspond to the upper half of the
  // circle centered at (1,1) with squared radius 3. We create the
  // circle with clockwise orientation, so the arc is directed from
  // (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
  Kernel::Point_2   c1 = Kernel::Point_2 (1, 1);
  Kernel::Circle_2  circ1 = Kernel::Circle_2 (c1, 3, CGAL::CLOCKWISE);
  CoordNT           one_minus_sqrt_3 = CoordNT (1, -1, 3);
  CoordNT           one_plus_sqrt_3 = CoordNT (1, 1, 3);
  Arc_point_2       s1 = Arc_point_2 (one_minus_sqrt_3, CoordNT (1));
  Arc_point_2       t1 = Arc_point_2 (one_plus_sqrt_3, CoordNT (1));
  curve1= Arc_section_2 (circ1, s1, t1);
  curve2= Arc_section_2 (circ1, s1, t1);
  //push the same semi circle again
  //curves.push_back (Arc_section_2 (circ1, s1, t1));

  //make x_monotone
  std::vector<CGAL::Object> X_monotone_curves;
  make_x_monotone_2(curve1, std::back_inserter(X_monotone_curves));
  make_x_monotone_2(curve2, std::back_inserter(X_monotone_curves));

  Arc_section_x_monotone_2 X_monotone_curve1, X_monotone_curve2, X_monotone_curve3 ;
  CGAL::assign( X_monotone_curve1, X_monotone_curves[ 0 ] );
  CGAL::assign( X_monotone_curve2, X_monotone_curves[ 1 ] );

  std::vector<CGAL::Object> Points_of_intersection;
  
  //intersect_2(X_monotone_curve1, X_monotone_curve2, std::back_inserter(Points_of_intersection));

 // Create a circular arc of the unit circle, directed clockwise from
  // (-1/2, sqrt(3)/2) to (1/2, sqrt(3)/2). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2  c6 = Kernel::Point_2 (0, 0);
  CoordNT          sqrt_3_div_2 = CoordNT (Number_type(0), Number_type(1,2), Number_type(3));
  Arc_point_2          s6 = Arc_point_2 (Number_type (-1, 2), sqrt_3_div_2);
  Arc_point_2          t6 = Arc_point_2 (Number_type (1, 2), sqrt_3_div_2);
  
  curve3 = Arc_section_2 (c6, 1, CGAL::CLOCKWISE, s6, t6);
  make_x_monotone_2(curve3, std::back_inserter(X_monotone_curves));
  CGAL::assign( X_monotone_curve2, X_monotone_curves[ 2 ] );

  Points_of_intersection.clear();
  //intersect_2(X_monotone_curve1, X_monotone_curve2, std::back_inserter(Points_of_intersection));
}

void check_compare_end_points_xy_2(Polycurve_arc_traits_2::Compare_endpoints_xy_2 compare_endpoints_xy_2, Polycurve_arc_traits_2::Make_x_monotone_2  make_x_monotone_2)
{
  Arc_section_2 curve1, curve2;

  // Create a circular arc that correspond to the upper half of the
  // circle centered at (1,1) with squared radius 3. We create the
  // circle with clockwise orientation, so the arc is directed from
  // (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
  Kernel::Point_2   c1 = Kernel::Point_2 (1, 1);
  Kernel::Circle_2  circ1 = Kernel::Circle_2 (c1, 3, CGAL::CLOCKWISE);
  CoordNT           one_minus_sqrt_3 = CoordNT (1, -1, 3);
  CoordNT           one_plus_sqrt_3 = CoordNT (1, 1, 3);
  Arc_point_2       s1 = Arc_point_2 (one_minus_sqrt_3, CoordNT (1));
  Arc_point_2       t1 = Arc_point_2 (one_plus_sqrt_3, CoordNT (1));
  curve1= Arc_section_2 (circ1, s1, t1);

  //make x_monotone
  std::vector<CGAL::Object> X_monotone_curves;
  make_x_monotone_2(curve1, std::back_inserter(X_monotone_curves));

  Arc_section_x_monotone_2 X_monotone_curve1, X_monotone_curve2 ;
  CGAL::assign( X_monotone_curve1, X_monotone_curves[ 0 ] );

  
  int res = compare_endpoints_xy_2(X_monotone_curve1);

  std::cout<< "The first result is: " << res << std::endl; 

  Kernel::Point_2   c2 = Kernel::Point_2 (1, 1);
  Kernel::Circle_2  circ2 = Kernel::Circle_2 (c2, 3, CGAL::COUNTERCLOCKWISE);
  Arc_point_2       t2 = Arc_point_2 (one_minus_sqrt_3, CoordNT (1));
  Arc_point_2       s2 = Arc_point_2 (one_plus_sqrt_3, CoordNT (1));
  curve2= Arc_section_2 (circ2, s2, t2);

  make_x_monotone_2(curve2, std::back_inserter(X_monotone_curves));
  CGAL::assign( X_monotone_curve2, X_monotone_curves[ 1 ] );

  res = compare_endpoints_xy_2(X_monotone_curve2);

  std::cout<< "The second result is: " << res << std::endl; 

}

void check_split( Polycurve_arc_traits_2::Split_2  split_2, Polycurve_arc_traits_2::Make_x_monotone_2  make_x_monotone_2)
{
  Arc_section_2 curve;

  // Create a circular arc that correspond to the upper half of the
  // circle centered at (1,1) with squared radius 3. We create the
  // circle with clockwise orientation, so the arc is directed from
  // (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
  Kernel::Point_2   c1 = Kernel::Point_2 (1, 1);
  Kernel::Circle_2  circ1 = Kernel::Circle_2 (c1, 3, CGAL::CLOCKWISE);
  CoordNT           one_minus_sqrt_3 = CoordNT (1, -1, 3);
  CoordNT           one_plus_sqrt_3 = CoordNT (1, 1, 3);
  Arc_point_2       s1 = Arc_point_2 (one_minus_sqrt_3, CoordNT (1));
  Arc_point_2       t1 = Arc_point_2 (one_plus_sqrt_3, CoordNT (1));
  curve= Arc_section_2 (circ1, s1, t1);

  //make x_monotone
  std::vector<CGAL::Object> X_monotone_curves;
  make_x_monotone_2(curve, std::back_inserter(X_monotone_curves));

  Arc_section_x_monotone_2 X_monotone_curve, split_x_monotone_curve1, split_x_monotone_curve2 ;
  CGAL::assign( X_monotone_curve, X_monotone_curves[ 0 ] );

  //split_2(X_monotone_curve, Kernel::Point_2::Kernel::Point_2(1, 4), split_x_monotone_curve1, split_x_monotone_curve2);

}

int main ()
{

  Polycurve_arc_traits_2 Polycurve_traits_2;

  //Compare_x_2
  //Polycurve_arc_traits_2::Compare_x_2 compare_x_2 = Polycurve_traits_2.compare_xy_2_object();

  //number of points
  //Polycurve_arc_traits_2::Number_of_points_2 num_of_points = Polycurve_traits_2.number_of_points_2_object();

  //construct min vertex
  Polycurve_arc_traits_2::Construct_min_vertex_2 cnst_min_vertex = Polycurve_traits_2.construct_min_vertex_2_object();

  //construct max vertex
  Polycurve_arc_traits_2::Construct_max_vertex_2 cnst_max_vertex_2 = Polycurve_traits_2.construct_max_vertex_2_object();

  //is vertical (return bool)
  Polycurve_arc_traits_2::Is_vertical_2 is_vertical = Polycurve_traits_2.is_vertical_2_object();

  //Compare y at x 2 (return comparison_result)
  Polycurve_arc_traits_2::Compare_y_at_x_2 cmp_y_at_x_2 = Polycurve_traits_2.compare_y_at_x_2_object();

  //compare y at x left
  Polycurve_arc_traits_2::Compare_y_at_x_left_2 cmp_y_at_x_left_2 = Polycurve_traits_2.compare_y_at_x_left_2_object();

  //compare y at x right
  Polycurve_arc_traits_2::Compare_y_at_x_right_2 cmp_y_at_x_right_2 = Polycurve_traits_2.compare_y_at_x_right_2_object();

  //equal_2
  Polycurve_arc_traits_2::Equal_2 equal_2 = Polycurve_traits_2.equal_2_object();

  //compare end points xy_2
  Polycurve_arc_traits_2::Compare_endpoints_xy_2 compare_endpoints_xy_2 = Polycurve_traits_2.compare_endpoints_xy_2_object();

  //construct opposite
  Polycurve_arc_traits_2::Construct_opposite_2  construct_opposite_2 = Polycurve_traits_2.construct_opposite_2_object();

  //make x_monotone
  Polycurve_arc_traits_2::Make_x_monotone_2  make_x_monotone_2 = Polycurve_traits_2.make_x_monotone_2_object();

  //push back
  Polycurve_arc_traits_2::Push_back_2  push_back_2 = Polycurve_traits_2.push_back_2_object();

  //push front
  Polycurve_arc_traits_2::Push_front_2  push_front_2 = Polycurve_traits_2.push_front_2_object();

  //split_2
  Polycurve_arc_traits_2::Split_2  split_2 = Polycurve_traits_2.split_2_object();

  //Intersect_2
  Polycurve_arc_traits_2::Intersect_2  intersect_2 = Polycurve_traits_2.intersect_2_object();

  //Are_mergable
  Polycurve_arc_traits_2::Are_mergeable_2  are_mergeable_2 = Polycurve_traits_2.are_mergeable_2_object();

  //Merge_2
  Polycurve_arc_traits_2::Merge_2  merge_2 = Polycurve_traits_2.merge_2_object();

  //construct_curve_2
  Polycurve_arc_traits_2::Construct_curve_2  construct_curve_2 = Polycurve_traits_2.construct_curve_2_object();

  //construct x_monotone curve_2
  Polycurve_arc_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = Polycurve_traits_2.construct_x_monotone_curve_2_object();



  check_equal(equal_2, make_x_monotone_2);
  check_intersect(make_x_monotone_2, intersect_2);
  check_compare_end_points_xy_2(compare_endpoints_xy_2, make_x_monotone_2);
  check_split(split_2, make_x_monotone_2);

  return 0;
}
#endif
