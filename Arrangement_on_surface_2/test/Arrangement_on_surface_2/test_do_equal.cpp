// Testing the do_equal function

#include <CGAL/Quotient.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_geometry_traits/Polycurve_2.h>
#include <CGAL/Arr_polyline_traits_2.h>

typedef CGAL::Quotient<int>                           Number_type;
typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2> Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;

int main ()
{
  bool are_equal;
  Traits_2 traits;
  Traits_2::Equal_2 equal = traits.equal_2_object();
  Traits_2::Construct_x_monotone_curve_2 ctr_xcv =
    traits.construct_x_monotone_curve_2_object();
  Point_2 points[5];
  X_monotone_curve_2 curve_1, curve_2;

  // Not equal case
  points[0] = Point_2(0,0);
  points[1] = Point_2(1,0);
  points[2] = Point_2(2,0);
  curve_1 = ctr_xcv(&points[0], &points[3]);
  points[0] = Point_2(0,0);
  points[1] = Point_2(1,0);
  points[2] = Point_2(3,0);
  curve_2 = ctr_xcv(&points[0], &points[3]);
  are_equal = equal(curve_1, curve_2);
  std::cout << "Answer1: " << are_equal<<std::endl;

  // Regular equal case
  points[0] = Point_2(0,0);
  points[1] = Point_2(1,0);
  points[2] = Point_2(2,0);
  curve_1 = ctr_xcv(&points[0], &points[3]);
  points[0] = Point_2(0,0);
  points[1] = Point_2(1,0);
  points[2] = Point_2(2,0);
  curve_2 = ctr_xcv(&points[0], &points[3]);
  are_equal = equal(curve_1, curve_2);
  std::cout << "Answer2: " << are_equal<<std::endl;

  // horizontal collinear case
  points[0] = Point_2(0,0);
  points[1] = Point_2(1,0);
  points[2] = Point_2(3,0);
  curve_1 = ctr_xcv(&points[0], &points[3]);
  points[0] = Point_2(0,0);
  points[1] = Point_2(1,0);
  points[2] = Point_2(2,0);
  points[3] = Point_2(3,0);
  curve_2 = ctr_xcv(&points[0], &points[4]);
  are_equal = equal(curve_1, curve_2);
  std::cout << "Answer3: " << are_equal<<std::endl;

  // collinear case
  points[0] = Point_2(0,0);
  points[1] = Point_2(2,2);
  points[2] = Point_2(4,4);
  curve_1 = ctr_xcv(&points[0], &points[3]);
  points[0] = Point_2(0,0);
  points[1] = Point_2(2,2);
  points[2] = Point_2(3,3);
  points[3] = Point_2(4,4);
  curve_2 = ctr_xcv(&points[0], &points[4]);
  are_equal = equal(curve_1, curve_2);
  std::cout << "Answer4: " << are_equal<<std::endl;

  // Vertical collinear case
  points[0] = Point_2(1,0);
  points[1] = Point_2(1,2);
  points[2] = Point_2(1,4);
  curve_1 = ctr_xcv(&points[0], &points[3]);
  points[0] = Point_2(1,0);
  points[1] = Point_2(1,2);
  points[2] = Point_2(1,3);
  points[3] = Point_2(1,4);
  curve_2 = ctr_xcv(&points[0], &points[4]);
  are_equal = equal(curve_1, curve_2);
  std::cout << "Answer5: " << are_equal<<std::endl;

  // Complicated not equal case
  points[0] = Point_2(0,0);
  points[1] = Point_2(1,1);
  points[2] = Point_2(2,2);
  points[3] = Point_2(4,2);
  curve_1 = ctr_xcv(&points[0], &points[4]);
  points[0] = Point_2(0,0);
  points[1] = Point_2(2,2);
  points[2] = Point_2(3,3);
  points[3] = Point_2(4,2);
  curve_2 = ctr_xcv(&points[0], &points[4]);
  are_equal = equal(curve_1, curve_2);
  std::cout << "Answer6: " << are_equal<<std::endl;

  // Complicated equal case
  points[0] = Point_2(0,0);
  points[1] = Point_2(1,1);
  points[2] = Point_2(2,2);
  points[3] = Point_2(4,2);
  curve_1 = ctr_xcv(&points[0], &points[4]);
  points[0] = Point_2(0,0);
  points[1] = Point_2(2,2);
  points[2] = Point_2(3,2);
  points[3] = Point_2(4,2);
  curve_2 = ctr_xcv(&points[0], &points[4]);
  are_equal = equal(curve_1, curve_2);
  std::cout << "Answer7: " << are_equal<<std::endl;

return 0;
}

