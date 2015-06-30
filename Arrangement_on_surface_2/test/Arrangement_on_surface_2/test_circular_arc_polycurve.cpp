
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

///////////////
//circle segment traits
//////////////
typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>     Arc_traits_2;
typedef CGAL::Arr_polyline_traits_2<Arc_traits_2>     Polycurve_arc_traits_2;


typedef Arc_traits_2::CoordNT                         CoordNT;
typedef Arc_traits_2::Point_2                         Point_2;
typedef Arc_traits_2::Curve_2                         Arc_section_2;
typedef Arc_traits_2::X_monotone_curve_2              Arc_x_monotone_section_2;

typedef Polycurve_arc_traits_2::X_monotone_curve_2    X_monotone_polycurve;
typedef Polycurve_arc_traits_2::Curve_2               Polycurve;
typedef Kernel::Circle_2                              Circle_2;

template<typename Curve>
bool check_compare_y_at_x_2(Curve& cv)
{
  Polycurve_arc_traits_2 traits;
  Polycurve_arc_traits_2::Compare_y_at_x_2 cmp_y_at_x_2 =
    traits.compare_y_at_x_2_object();

  //create points
  Polycurve_arc_traits_2::Point_2 point_above_line =
    Polycurve_arc_traits_2::Point_2(4,4),
    point_below_line = Polycurve_arc_traits_2::Point_2(3,1),
    point_on_line = Polycurve_arc_traits_2::Point_2(0,1);

  CGAL::Comparison_result result;

  result =  cmp_y_at_x_2(point_above_line, cv);
  std::cout << "Compare_y_at_x_2:: for point above the curve computed Answer is:  "
            << (result == CGAL::SMALLER ? "Below":
               (result == CGAL::LARGER ? "Above" : "On-line")) << std::endl;

  result =  cmp_y_at_x_2(point_below_line, cv);
  std::cout << "Compare_y_at_x_2:: for point below the curve computed Answer is:  "
            << (result == CGAL::SMALLER ? "Below":
               (result == CGAL::LARGER ? "Above" : "On-line")) << std::endl;

  result =  cmp_y_at_x_2(point_on_line, cv);
  std::cout << "Compare_y_at_x_2:: for point on the curve computed Answer is:  "
            << (result == CGAL::SMALLER ? "Below":
               (result == CGAL::LARGER ? "Above" : "On-line")) << std::endl;

  return true;
}

template<typename Curve>
void check_intersect(const Curve& cv1, const Curve& cv2)
{
  // Polycurve_arc_traits_2 traits;

  // std::vector<CGAL::Object> object_vec;
  // traits.intersect_2_object()(cv1, cv2, object_vec);
  // std::cout<< "number of intersections is: " << object_vec.size();
}

template<typename Curve>
void check_make_x_monotone(Curve cv)
{
   Polycurve_arc_traits_2 traits;
   std::vector<CGAL::Object> object_vec;

   traits.make_x_monotone_2_object()(cv, std::back_inserter(object_vec));
   std::cout << "Number of x-monotone curves: "
             << object_vec.size() << std::endl;
}

template<typename Curve>
void check_trim(Curve& xcv, int sx, int sy, int tx, int ty)
{
  Polycurve_arc_traits_2 traits;
  Point_2 source(sx, sy);
  Point_2 target(tx, ty);

  Polycurve_arc_traits_2::Trim_2  trim_polycurve = traits.trim_2_object();
  X_monotone_polycurve trimmed_curve = trim_polycurve(xcv, source, target);

  std::cout << "polycurvecurve: " << xcv << std::endl<<std::endl;
  std::cout << "Trimmed curve: " << trimmed_curve << std::endl;
}

int main(int argc, char* argv[])
{
  Polycurve_arc_traits_2 traits;

  std::vector<Arc_section_2> curves;
  std::vector<Arc_x_monotone_section_2> x_curves;

  // Create a circular arc of the circle, directed clockwise from
  // (-1, 0) to (1, 0). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2 c6 = Kernel::Point_2(0, 0);
  Point_2 s6 = Point_2(Number_type(-1, 1), Number_type(0, 1));
  Point_2 t6 = Point_2(Number_type(1, 1), Number_type(0, 1));
  Arc_section_2 circ_arc1(c6, 1, CGAL::CLOCKWISE, s6, t6);

  curves.push_back(circ_arc1);

  // Create a circular arc of the unit circle, directed clockwise from
  // (1, 0) to (3, 0). Note that we orient the
  // supporting circle accordingly.
  Kernel::Point_2 c1 = Kernel::Point_2(3, 0);
  Point_2 s1 = Point_2(Number_type(1, 1), Number_type(0, 1));
  Point_2 t1 = Point_2(Number_type(5, 1), Number_type(0, 1));
  Arc_section_2 circ_arc2(c1, 2, CGAL::CLOCKWISE, s1, t1);
  curves.push_back(circ_arc2);

  Circle_2 circ = Circle_2(c6, 1, CGAL::CLOCKWISE);
  Arc_x_monotone_section_2 xc1(circ, s6, t6, CGAL::CLOCKWISE);
  x_curves.push_back(xc1);

  Circle_2 circ1 = Circle_2(c1, 4, CGAL::CLOCKWISE);
  Arc_x_monotone_section_2 xc2(circ1, s1, t1, CGAL::CLOCKWISE);
  x_curves.push_back(xc2);

  // Polycurve circ_arc_polycurve =
  //   traits.construct_curve_2_object()(curves.begin(), curves.end());
  X_monotone_polycurve x_polycurve_1 =
    traits.construct_x_monotone_curve_2_object()(x_curves.begin(),
                                                 x_curves.end());

  Kernel::Point_2 cen = Kernel::Point_2(-2, 0);
  Circle_2 circ2 = Circle_2(cen, 2, CGAL::CLOCKWISE);
  Point_2 s2 = Point_2(Number_type(-4, 1), Number_type(0, 1));
  Point_2 t2 = Point_2(Number_type(0, 1), Number_type(0, 1));
  Arc_x_monotone_section_2 xc3(circ2, s2, t2, CGAL::CLOCKWISE);
  x_curves.clear();
  x_curves.push_back(xc3);
  X_monotone_polycurve x_polycurve_2 =
    traits.construct_x_monotone_curve_2_object()(x_curves.begin(),
                                                 x_curves.end());

  //testing for arc construction from two points.
  //Arc_x_monotone_section_2 x_segment(Kernel::Point_2(0, 0),
  //                                   Kernel::Point_2(2, 0));
  //x_curves.clear();
  //x_curves.push_back(x_segment);
  //X_monotone_polycurve x_polycurve_3 =
  //traits.construct_x_monotone_curve_2_object()(x_curves.begin(),
  //                                             x_curves.end());
  //std::cout<< "x_polycurve_3: " << x_polycurve_3 << std::endl;

  //Another polycurve
  curves.clear();

  Kernel::Point_2  center = Kernel::Point_2(-10, 10);
  Point_2 source = Point_2(Number_type(-10, 1), Number_type(13, 1));
  Point_2 target = Point_2(Number_type(-7, 1), Number_type(10, 1));
  Arc_section_2 circ_arc(center, 3, CGAL::CLOCKWISE, source, target);
  curves.push_back(circ_arc);

  center = Kernel::Point_2(-20, 10);
  source = Point_2(Number_type(-7, 1), Number_type(10, 1));
  target = Point_2(Number_type(-20, 1), Number_type(23, 1));
  circ_arc =  Arc_section_2(center, 13, CGAL::CLOCKWISE, source, target);
  curves.push_back(circ_arc);

  center = Kernel::Point_2(-20, 25);
  source = Point_2(Number_type(-20, 1), Number_type(23, 1));
  target = Point_2(Number_type(-20, 1), Number_type(27, 1));
  circ_arc =  Arc_section_2(center, 2, CGAL::CLOCKWISE, source, target);
  curves.push_back(circ_arc);

  Polycurve curve_1 =
    traits.construct_curve_2_object()(curves.begin(), curves.end());

  //////////////////////
  //Functor testing
  //////////////////////

  // check_compare_y_at_x_2(x_polycurve_1);
  // check_intersect(x_polycurve_1, x_polycurve_2);
  //check_make_x_monotone(curve_1);

  //checking if the cgal_assertion for curve construction for two points work
  //or not.
  //Point_2 push_back_point(Number_type(10, 1), Number_type(0, 1));
  //traits.push_back_2_object()(x_polycurve_1, push_back_point);

  //  //checking for trim.
  // Arc_traits_2 arc_traits;
  // source = Point_2(Number_type(1, 1), Number_type(0, 1));
  // target = Point_2(Number_type(3, 1), Number_type(2, 1));
  // // source = Point_2(Number_type(2, 1), Number_type(-2, 1));
  // // target = Point_2(Number_type(3, 1), Number_type(4, 1));
  // std::cout << " curve is : " << xc2 << std::endl;
  // Arc_x_monotone_section_2 trimmed_curve =
  //   arc_traits.trim_2_object()(xc2, source, target);
  // std::cout << "trimmed conic curve is : " << trimmed_curve << std::endl;

  check_trim(x_polycurve_1, atoi(argv[1]), atoi(argv[2]),
             atoi(argv[3]), atoi(argv[4]));
  std::cout << std::endl;

  return 0;
}
#endif
