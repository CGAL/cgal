// Testing the do_equal function

#include <CGAL/basic.h>
#ifndef CGAL_USE_CORE
#include <iostream>
int main()
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
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <boost/type_traits/is_same.hpp>

////////////////////
//conic traits
////////////////////
typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef Rat_kernel::Point_2                           Rat_point_2;
typedef Rat_kernel::Segment_2                         Rat_segment_2;
typedef Rat_kernel::Circle_2                          Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                      Conic_traits_2;
typedef Conic_traits_2::Point_2                       Conic_point_2;
typedef Conic_traits_2::Curve_2                       Conic_curve_2;
typedef Conic_traits_2::X_monotone_curve_2            Conic_x_monotone_curve_2;
typedef CGAL::Arr_polyline_traits_2<Conic_traits_2>   Polycurve_conic_traits_2;
typedef Polycurve_conic_traits_2::X_monotone_curve_2  Pc_x_monotone_curve_2;
// typedef Polycurve_conic_traits_2::Point_2          polypoint;

// typedef CGAL::Arr_polyline_traits_2<
//                CGAL::Arr_conic_traits_2<CGAL::Cartesian<BigRat>,
//                CGAL::Cartesian<Expr>,
//                CGAL::CORE_algebraic_number_traits>
//              >::Point_2   test_point_2;

void check_equal()
{
  bool are_equal;

  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Equal_2 equal = traits.equal_2_object();
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();

  //create some curves
  Conic_point_2 ps1(Rational(1,4), 4);
  Conic_point_2 pt1(2, Rational(1,2));
  Conic_curve_2 c1(0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE, ps1, pt1);

  Conic_point_2 ps2(Rational(1,4), 4);
  Conic_point_2 pt2(2, Rational(1,2));
  Conic_curve_2 c2(0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE, ps2, pt2);

  Rat_point_2 ps3(Rational(1,4), 4);
  Rat_point_2 pmid3(Rational(3,2), 2);
  Rat_point_2 pt3(2, Rational(1,3));
  Conic_curve_2 c3(ps3, pmid3, pt3);

  Rat_point_2 ps4(1, 5);
  Rat_point_2 pmid4(Rational(3,2), 3);
  Rat_point_2 pt4(3, Rational(1,3));
  Conic_curve_2 c4(ps4, pmid4, pt4);

  // //make x_monotone
  Polycurve_conic_traits_2::X_monotone_curve_2 xmc1 =
    construct_x_monotone_curve_2(c1);
  Polycurve_conic_traits_2::X_monotone_curve_2 xmc2 =
    construct_x_monotone_curve_2(c2);
  Polycurve_conic_traits_2::X_monotone_curve_2 xmc3 =
    construct_x_monotone_curve_2(c3);
  Polycurve_conic_traits_2::X_monotone_curve_2 xmc4 =
    construct_x_monotone_curve_2(c4);

  are_equal = equal(xmc1, xmc2);
  std::cout << "Two equal conic arcs are computed as:  "
            << ((are_equal) ? "equal" : "Not equal") << std::endl;

  are_equal = equal(xmc3, xmc2);
  std::cout << "Two un-equal conic arcs are computed as:  "
            << ((are_equal) ? "equal" : "Not equal") << std::endl;

  are_equal = equal(xmc3, xmc4);
  std::cout << "Two un-equal conic arcs are computed as:  "
            << ((are_equal) ? "equal" : "Not equal") << std::endl;
 }

 template <typename curve_type>
 void check_intersect(curve_type &xcv1, curve_type &xcv2)
 {
  Polycurve_conic_traits_2 traits;
  std::vector<CGAL::Object> intersection_points;
  traits.intersect_2_object()(xcv1, xcv2, std::back_inserter(intersection_points));
  std::cout<< "Number of intersection Points: " << intersection_points.size()
           << std::endl;

  //dynamic cast the cgal_objects
  // std::vector< std::pair<Polycurve_conic_traits_2::Point_2,
  //                        Polycurve_conic_traits_2::Multiplicity> > pm_vector;
  // for(int i=0; i<intersection_points.size(); i++)
  // {
  //   std::pair<Polycurve_conic_traits_2::Point_2,
  //             Polycurve_conic_traits_2::Multiplicity> pm =
  //   CGAL::object_cast<std::pair<Polycurve_conic_traits_2::Point_2,
  //                               Polycurve_conic_traits_2::Multiplicity> >
  //                              (&(intersection_points[i]));
  //   pm_vector.push_back(pm);
  // }
 }

void check_compare_end_points_xy_2()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Compare_endpoints_xy_2 compare_endpoints_xy_2 =
    traits.compare_endpoints_xy_2_object();

  //create some curves
  Conic_point_2 ps1(Rational(1,4), 4);
  Conic_point_2 pt1(2, Rational(1,2));
  Conic_curve_2 c1(0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE, ps1, pt1);

  // Insert a parabolic arc that is supported by a parabola y = -x^2
  // (or: x^2 + y = 0) and whose endpoints are (-sqrt(3), -3) ~ (-1.73, -3)
  // and (sqrt(2), -2) ~ (1.41, -2). Notice that since the x-coordinates
  // of the endpoints cannot be acccurately represented, we specify them
  // as the intersections of the parabola with the lines y = -3 and y = -2.
  // Note that the arc is clockwise oriented.
  Conic_curve_2
    c2 = Conic_curve_2(1, 0, 0, 0, 1, 0,         // The parabola.
                       CGAL::CLOCKWISE,
                       Conic_point_2(-1.73, -3), // Approximation of the source.
                       0, 0, 0, 0, 1, 3,         // The line: y = -3.
                       Conic_point_2(1.41, -2),  // Approximation of the target.
                       0, 0, 0, 0, 1, 2);        // The line: y = -2.
  CGAL_assertion(c2.is_valid());

  //make polyline x-monotone curves
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 =
    construct_x_monotone_curve_2(c1);
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 =
    construct_x_monotone_curve_2(c2);

  CGAL::Comparison_result res = compare_endpoints_xy_2(polyline_xmc1);
  std::cout << "compare_end_points_xy_2 for counterclockwise curve: "
            << (res == CGAL::SMALLER ? "SMALLER":
               (res == CGAL::LARGER ? "LARGER" : "EQUAL")) << std::endl;

  res = compare_endpoints_xy_2(polyline_xmc2);
  std::cout<< "compare_end_points_xy_2 for clockwise curve: "
           << (res == CGAL::SMALLER ? "SMALLER":
               (res == CGAL::LARGER ? "LARGER" : "EQUAL")) << std::endl;
}

template <typename Curve_type>
void check_split(Curve_type &xcv1, Curve_type &xcv2)
{
  Polycurve_conic_traits_2 traits;

  //split x poly-curves

  Conic_curve_2 c6(1,1,0,6,-26,162,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(-7), Algebraic(13)),
                   Conic_point_2(Algebraic(-3), Algebraic(9)));
  Conic_curve_2 c7(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(-3), Algebraic(9)),
                   Conic_point_2(Algebraic(0), Algebraic(0)));
  Conic_curve_2 c8(0,1,0,-1,0,0, CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(0), Algebraic(0)),
                   Conic_point_2(Algebraic(4), Algebraic(-2)));

  Conic_x_monotone_curve_2 xc6(c6);
  Conic_x_monotone_curve_2 xc7(c7);
  Conic_x_monotone_curve_2 xc8(c8);
  std::vector<Conic_x_monotone_curve_2> xmono_conic_curves_2;

  xmono_conic_curves_2.push_back(xc6);
  xmono_conic_curves_2.push_back(xc7);
  Pc_x_monotone_curve_2 split_expected_1 =
    traits.construct_x_monotone_curve_2_object()(xmono_conic_curves_2.begin(),
                                                 xmono_conic_curves_2.end());

  xmono_conic_curves_2.clear();
  xmono_conic_curves_2.push_back(xc8);
  Pc_x_monotone_curve_2 split_expected_2 =
    traits.construct_x_monotone_curve_2_object()(xmono_conic_curves_2.begin(),
                                                 xmono_conic_curves_2.end());


  Polycurve_conic_traits_2::X_monotone_curve_2 split_curve_1, split_curve_2;
  Polycurve_conic_traits_2::Point_2
    point_of_split = Polycurve_conic_traits_2::Point_2(0,0);

  //Split functor
  traits.split_2_object()(xcv2, point_of_split, split_curve_1, split_curve_2);

  bool split_1_chk = traits.equal_2_object()(split_curve_1, split_expected_1);
  bool split_2_chk = traits.equal_2_object()(split_curve_2, split_expected_2);

  if(split_1_chk && split_2_chk)
    std::cout << "Split is working fine" << std::endl;
  else
    std::cout << "Something is wrong with split" << std::endl;
}

void check_is_vertical()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Is_vertical_2 is_vertical =
    traits.is_vertical_2_object();

   //create a curve
  Rat_point_2 ps1(1, 10);
  Rat_point_2 pmid1(5, 4);
  Rat_point_2 pt1(10, 1);
  Conic_curve_2 c1(ps1, pmid1, pt1);

  //make x-monotone curve
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 =
    construct_x_monotone_curve_2(c1);

  bool result = is_vertical(polyline_xmc1);
  std::cout << "Is_verticle:: Expected first result is not vertivle: Computed: "
            << ((result)? "vertical" : "not vertical") << std::endl;
}

/*! */
template <typename stream>
bool read_orientation(stream& is, CGAL::Orientation& orient)
{
  int i_orient;
  is >> i_orient;
  orient = (i_orient > 0) ? CGAL::COUNTERCLOCKWISE :
    (i_orient < 0) ? CGAL::CLOCKWISE : CGAL::COLLINEAR;
  return true;
}

/*! */
template <typename stream>
bool read_app_point(stream& is, Conic_point_2& p)
{
  //waqar: original
  double x, y;
  is >> x >> y;
  p = Conic_point_2(Algebraic(x), Algebraic(y));
  return true;

  //waqar: my modification
  // long int rat_x_num, rat_x_den, rat_y_num, rat_y_den;
  // is >> rat_x_num >> rat_x_den >> rat_y_num >> rat_y_den;

  // Basic_number_type x(rat_x), y(rat_y);
  // p = Conic_point_2(Rational(rat_x_num, rat_x_den),
  //                   Rational(rat_y_num, rat_y_den));
  // return true;
}

/*! */
template <typename stream>
bool read_orientation_and_end_points(stream& is, CGAL::Orientation& orient,
                                     Conic_point_2& source,
                                     Conic_point_2& target)
{
  // Read the orientation.
  if (!read_orientation(is, orient)) return false;

  // Read the end points of the arc and create it.
  if (!read_app_point(is, source)) return false;
  if (!read_app_point(is, target)) return false;
  return true;
}

/*! */
template <typename stream, typename Curve>
bool read_general_arc(stream& is, Curve& cv)
{
  // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
  Rational r, s, t, u, v, w;                // The conic coefficients.
  is >> r >> s >> t >> u >> v >> w;
    // Read the orientation.
  int i_orient = 0;
  is >> i_orient;
  CGAL::Orientation orient = (i_orient > 0) ? CGAL::COUNTERCLOCKWISE :
    (i_orient < 0) ? CGAL::CLOCKWISE : CGAL::COLLINEAR;

  // Read the approximated source, along with a general conic
  // <r_1,s_1,t_1,u_1,v_1,w_1> whose intersection with <r,s,t,u,v,w>
  // defines the source.
  Conic_point_2 app_source;
  if (!read_app_point(is, app_source)) return false;
  Rational r1, s1, t1, u1, v1, w1;
  is >> r1 >> s1 >> t1 >> u1 >> v1 >> w1;

  // Read the approximated target, along with a general conic
  // <r_2,s_2,t_2,u_2,v_2,w_2> whose intersection with <r,s,t,u,v,w>
  // defines the target.
  Conic_point_2 app_target;
  if (!read_app_point(is, app_target)) return false;

  Rational r2, s2, t2, u2, v2, w2;
  is >> r2 >> s2 >> t2 >> u2 >> v2 >> w2;

  std::cout << "line is: " << r << s << t << u << v << w << i_orient
            << r1 << s1 << t1 << u1 << v1 << w1
            << r2 << s2 << t2 << u2 << v2 << w2 << std::endl;

  // Create the conic arc.
  cv = Curve(r, s, t, u, v, w, orient,
             app_source, r1, s1, t1, u1, v1, w1,
             app_target, r2, s2, t2, u2, v2, w2);
  return true;
}

/*! */
template <typename stream, typename Curve>
bool read_general_conic(stream& is, Curve& cv)
{
  // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
  Rational r, s, t, u, v, w;
  is >> r >> s >> t >> u >> v >> w;
  // Create a full conic (should work only for ellipses).
  cv = Curve(r, s, t, u, v, w);
  return true;
}

// /*! */
template <typename stream, typename Curve>
bool read_general_curve(stream& is, Curve& cv)
{
  Rational r, s, t, u, v, w;                // The conic coefficients.
  // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
  is >> r >> s >> t >> u >> v >> w;
  CGAL::Orientation orient;
  Conic_point_2 source, target;
  if (!read_orientation_and_end_points(is, orient, source, target))
    return false;

  // Create the conic (or circular) arc.
  // std::cout << "arc coefficients: " << r << " " << s << " " << t << " "
  //           << u << " " << v << " " << w << std::endl;
  // std::cout << "Read Points : " << source.x() << " " << source.y() << " "
  //           << target.x() << " " << target.y() << std::endl;
  cv = Curve(r, s, t, u, v, w, orient, source, target);
  return true;
}
std::istream& skip_comments(std::istream& is, std::string& line)
{
  while (std::getline(is, line))
    if (!line.empty() && (line[0] != '#')) break;
  return is;
}

bool check_compare_y_at_x_2()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Compare_y_at_x_2 cmp_y_at_x_2 =
    traits.compare_y_at_x_2_object();
  //polycurve constructors
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_mono_polycurve = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Construct_curve_2  construct_polycurve =
    traits.construct_curve_2_object();

   //create a curve
  Rat_point_2 ps1(1, 10);
  Rat_point_2 pmid1(5, 4);
  Rat_point_2 pt1(10, 1);
  Conic_curve_2 c1(ps1, pmid1, pt1);

   //create a curve
  Rat_point_2 ps2(10, 1);
  Rat_point_2 pmid2(15, 5);
  Rat_point_2 pt2(20, 10);
  Conic_curve_2 c2(ps2, pmid2, pt2);

  Conic_curve_2 c3(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(0), Algebraic(0)),
                   Conic_point_2(Algebraic(3), Algebraic(9)));
  Conic_curve_2 c4(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(3), Algebraic(9)),
                   Conic_point_2(Algebraic(5), Algebraic(25)));

  std::vector<Conic_curve_2> conic_curves, conic_curves_2, Conic_curves_3;
  conic_curves.push_back(c1);
  conic_curves.push_back(c2);

  //conic_curves_2.push_back(c3);
  //conic_curves_2.push_back(c4);

  Conic_x_monotone_curve_2 xc1(c1);
  Conic_x_monotone_curve_2 xc2(c2);
  Conic_x_monotone_curve_2 xc3(c3);
  Conic_x_monotone_curve_2 xc4(c4);

  std::vector<Conic_x_monotone_curve_2> xmono_conic_curves, xmono_conic_curves_2;
  /* VERY IMPORTANT
   * For efficiency reasons, we recommend users not to construct x-monotone
   * conic arc directly, but rather use the Make_x_monotone_2 functor supplied
   * by the conic-arc traits class to convert conic curves to x-monotone curves.
   */
  xmono_conic_curves.push_back(xc1);
  xmono_conic_curves.push_back(xc2);
  xmono_conic_curves_2.push_back(xc3);
  xmono_conic_curves_2.push_back(xc4);

  //construct x-monotone poly-curve
  Polycurve_conic_traits_2::X_monotone_curve_2 conic_x_mono_polycurve =
    construct_x_mono_polycurve(xmono_conic_curves.begin(),
                               xmono_conic_curves.end());
  Polycurve_conic_traits_2::X_monotone_curve_2 conic_x_mono_polycurve_2 =
    construct_x_mono_polycurve(xmono_conic_curves_2.begin(),
                               xmono_conic_curves_2.end());

  //construct poly-curve
  Polycurve_conic_traits_2::Curve_2 conic_polycurve =
    construct_polycurve(conic_curves.begin(), conic_curves.end());
  //Polycurve_conic_traits_2::Curve_2 conic_polycurve_2 =
  //  construct_polycurve(conic_curves_2.begin(), conic_curves_2.end());

  //make x-monotone curve
  //Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 =
  //  construct_x_monotone_curve_2(c1);

  //create points
  Polycurve_conic_traits_2::Point_2
    point_above_line = Polycurve_conic_traits_2::Point_2(2,10),
    point_below_line = Polycurve_conic_traits_2::Point_2(4,7),
    point_on_line = Polycurve_conic_traits_2::Point_2(2,4);

  CGAL::Comparison_result result;

  result =  cmp_y_at_x_2(point_above_line, conic_x_mono_polycurve_2);
  std::cout << "Compare_y_at_x_2:: for point above the curve computed Answer is:  "
            << (result == CGAL::SMALLER ? "Below":
               (result == CGAL::LARGER ? "Above" : "On-line")) << std::endl;

  result =  cmp_y_at_x_2(point_below_line, conic_x_mono_polycurve_2);
  std::cout << "Compare_y_at_x_2:: for point below the curve computed Answer is:  "
            << (result == CGAL::SMALLER ? "Below":
               (result == CGAL::LARGER ? "Above" : "On-line")) << std::endl;

  result =  cmp_y_at_x_2(point_on_line, conic_x_mono_polycurve_2);
  std::cout << "Compare_y_at_x_2:: for point on the curve computed Answer is:  "
            << (result == CGAL::SMALLER ? "Below":
               (result == CGAL::LARGER ? "Above" : "On-line")) << std::endl;

  return true;
}

void check_are_mergable()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Are_mergeable_2  are_mergeable_2 =
    traits.are_mergeable_2_object();

   //create a curve
  Rat_point_2 ps1(1, 10);
  Rat_point_2 pmid1(5, 4);
  Rat_point_2 pt1(10, 1);
  Conic_curve_2 c1(ps1, pmid1, pt1);

  Rat_point_2 ps2(10, 1);
  Rat_point_2 pmid2(15, 14);
  Rat_point_2 pt2(20, 20);
  Conic_curve_2 c2(ps2, pmid2, pt2);

  Rat_point_2 ps3(Rational(1,4), 4);
  Rat_point_2 pmid3(Rational(3,2), 2);
  Rat_point_2 pt3(2, Rational(1,3));
  Conic_curve_2 c3(ps3, pmid3, pt3);

  //construct x-monotone curve(compatible with polyline class)
   Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 =
     construct_x_monotone_curve_2(c1);
   Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 =
     construct_x_monotone_curve_2(c2);
   Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc3 =
     construct_x_monotone_curve_2(c3);

  bool result = are_mergeable_2(polyline_xmc1, polyline_xmc2);
  std::cout << "Are_mergable:: Mergable x-monotone polycurves are Computed as: "
            << ((result)? "Mergable" : "Not-Mergable") << std::endl;

  result = are_mergeable_2(polyline_xmc1, polyline_xmc3);
  std::cout << "Are_mergable:: Non-Mergable x-monotone polycurves are Computed as: "
            << ((result)? "Mergable" : "Not-Mergable") << std::endl;
}

void check_merge_2()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Merge_2  merge_2 = traits.merge_2_object();

  //create a curve
  Rat_point_2 ps1(1, 10);
  Rat_point_2 pmid1(5, 4);
  Rat_point_2 pt1(10, 1);
  Conic_curve_2 c1(ps1, pmid1, pt1);

  Rat_point_2 ps2(10, 1);
  Rat_point_2 pmid2(15, 14);
  Rat_point_2 pt2(20, 20);
  Conic_curve_2 c2(ps2, pmid2, pt2);

//construct x-monotone curve (compatible with polyline class)
 Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 =
   construct_x_monotone_curve_2(c1);
 Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 =
   construct_x_monotone_curve_2(c2);

 Polycurve_conic_traits_2::X_monotone_curve_2 merged_xmc;

 merge_2(polyline_xmc1, polyline_xmc2, merged_xmc);
 std::cout<< "Merge_2:: Mergable x-monotone curves merged successfully"
          << std:: endl;
}

void check_construct_opposite()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Construct_opposite_2 construct_opposite_2 =
    traits.construct_opposite_2_object();

  //create a curve
  Rat_point_2 ps1(1, 10);
  Rat_point_2 pmid1(5, 4);
  Rat_point_2 pt1(10, 1);
  Conic_curve_2     c1(ps1, pmid1, pt1);

  //construct x-monotone curve (compatible with polyline class)
 Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 =
   construct_x_monotone_curve_2(c1);
 Polycurve_conic_traits_2::X_monotone_curve_2 polyline_opposite_curve =
   construct_opposite_2(polyline_xmc1);

 std::cout<< "Construct_opposite_2:: Opposite curve created";
}

void check_compare_y_at_x_right()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Compare_y_at_x_right_2 cmp_y_at_x_right_2 =
    traits.compare_y_at_x_right_2_object();

  //create constructing curves
  Rat_point_2 ps2(1, 10);
  Rat_point_2 pmid2(5, 4);
  Rat_point_2 pt2(10, 1);
  Conic_curve_2 c1(ps2, pmid2, pt2);

  Rat_point_2 ps3(10, 1);
  Rat_point_2 pmid3(5, 4);
  Rat_point_2 pt3(1, 10);
  Conic_curve_2 c2(ps3, pmid3, pt3);

  //construct x-monotone curve (compatible with polyline class)
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 =
    construct_x_monotone_curve_2(c1);
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 =
    construct_x_monotone_curve_2(c2);
  Polycurve_conic_traits_2::Point_2 intersection_point =
    Polycurve_conic_traits_2::Point_2(5,4);

  CGAL::Comparison_result result;
  result = cmp_y_at_x_right_2(polyline_xmc1, polyline_xmc2, intersection_point);
  std::cout << "Compare_y_at_x_right:: Expected Answer: equal, Computed answer:  "
            << (result == CGAL::SMALLER ? "smaller":
               (result == CGAL::LARGER ? "Larger" : "equal")) << std::endl;
}

void check_compare_y_at_x_left()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Compare_y_at_x_left_2 cmp_y_at_x_left_2 =
    traits.compare_y_at_x_left_2_object();

  //create constructing curves
  Rat_point_2 ps2(1, 10);
  Rat_point_2 pmid2(5, 4);
  Rat_point_2 pt2(10, 1);
  Conic_curve_2 c1(ps2, pmid2, pt2);

  Rat_point_2 ps3(10, 1);
  Rat_point_2 pmid3(5, 4);
  Rat_point_2 pt3(1, 10);
  Conic_curve_2 c2(ps3, pmid3, pt3);

  //construct x-monotone curve(compatible with polyline class)
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 =
    construct_x_monotone_curve_2(c1);
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 =
    construct_x_monotone_curve_2(c2);
  Polycurve_conic_traits_2::Point_2 intersection_point =
    Polycurve_conic_traits_2::Point_2(5,4);

  CGAL::Comparison_result result;

  result = cmp_y_at_x_left_2(polyline_xmc1, polyline_xmc2, intersection_point);
  std::cout << "Compare_y_at_x_left:: Expected Answer: equal, Computed answer:  "
            << (result == CGAL::SMALLER ? "smaller":
               (result == CGAL::LARGER ? "Larger" : "equal")) << std::endl;
}

template <typename Curve_type>
void check_make_x_monotne_curve(Curve_type c1)
{
  Polycurve_conic_traits_2 traits;
  std::vector<CGAL::Object> obj_vec;
  traits.make_x_monotone_2_object()(c1, std::back_inserter(obj_vec));

  std::cout << "The polycurve is: " << c1 << std::endl;

  std::cout<< "The poly curve have been split into " << obj_vec.size()
           << " polycurves" << std::endl;

  //const Pc_x_monotone_curve_2 *split_curve_1 =
  //  CGAL::object_cast<Pc_x_monotone_curve_2> (&(obj_vec[0]));
  //const Pc_x_monotone_curve_2 *split_curve_2 =
  //  CGAL::object_cast<Pc_x_monotone_curve_2> (&(obj_vec[1]));

  //std::cout << "The split curve 1 is: " << *split_curve_1 << std::endl;
  //std::cout << "The split curve 2 is: " << *split_curve_2 << std::endl;
}

template<typename Curve, typename Segment>
void check_push_front(Curve base_curve, Segment curve_tobe_pushed)
{
  Polycurve_conic_traits_2 traits;

  std::cout << "Base curve: " << base_curve << std::endl;
  std::cout << "push curve: " << curve_tobe_pushed << std::endl;

  traits.push_front_2_object()(base_curve, curve_tobe_pushed);
  std::cout << "result curve: " << base_curve << std::endl;
}

template<typename Curve, typename Segment>
void check_push_back(Curve& base_curve, Segment curve_tobe_pushed)
{
  Polycurve_conic_traits_2 traits;

  std::cout << "Base curve: " << base_curve << std::endl;
  std::cout << "push curve: " << curve_tobe_pushed << std::endl;

  traits.push_back_2_object()(base_curve, curve_tobe_pushed);

  std::cout << "result curve: " << base_curve << std::endl;
}

template<typename Segment>
void check_compare_x_2(const Segment& seg1, const Segment& seg2)
{
  Polycurve_conic_traits_2 traits;
  CGAL::Comparison_result result;

  result = traits.compare_x_2_object()(seg1, CGAL::ARR_MIN_END, seg2,
                                       CGAL::ARR_MIN_END);
  std::cout << "Compare_x_2:: Expected Answer: Larger, Computed answer:  "
            << (result == CGAL::SMALLER ? "smaller":
               (result == CGAL::LARGER ? "Larger" : "equal")) << std::endl;

  result = traits.compare_x_2_object()(seg1, CGAL::ARR_MIN_END, seg2,
                                       CGAL::ARR_MAX_END);
  std::cout << "Compare_x_2:: Expected Answer: Equal, Computed answer:  "
            << (result == CGAL::SMALLER ? "smaller":
               (result == CGAL::LARGER ? "Larger" : "equal")) << std::endl;
}

template<typename Curve>
void check_compare_points(Curve& cv)
{
  Polycurve_conic_traits_2 traits;
  CGAL::Arr_parameter_space result =
    traits.parameter_space_in_x_2_object()(cv, CGAL::ARR_MAX_END);
}

template <typename curve>
void check_trim(curve& xcv, int sx, int sy, int tx, int ty)
{
  Polycurve_conic_traits_2 traits;

  // Conic_point_2 source(Algebraic(-16), Algebraic(-4));
  // Conic_point_2 target(Algebraic(4), Algebraic(16));
  Conic_point_2 source(sx, sy);
  Conic_point_2 target(tx, ty);

  Polycurve_conic_traits_2::Trim_2  trim_polycurve = traits.trim_2_object();
  Pc_x_monotone_curve_2 trimmed_curve = trim_polycurve(xcv, source, target);

  std::cout << "polycurvecurve: " << xcv << std::endl<<std::endl;
  std::cout << "Trimmed curve: " << trimmed_curve << std::endl;

}

int main(int argc, char* argv[])
{
  Polycurve_conic_traits_2 traits;
    //polycurve constructors
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2
    construct_x_mono_polycurve = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Construct_curve_2  construct_polycurve =
    traits.construct_curve_2_object();

   //create a curve

  Conic_curve_2 c3(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(0), Algebraic(0)),
                   Conic_point_2(Algebraic(3), Algebraic(9)));
  Conic_curve_2 c4(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(3), Algebraic(9)),
                   Conic_point_2(Algebraic(5), Algebraic(25)));
  Conic_curve_2 c5(0,1,0,1,0,0, CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(-25), Algebraic(-5)),
                   Conic_point_2(Algebraic(0), Algebraic(0)));

  Conic_curve_2 c6(1,1,0,6,-26,162,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(-7), Algebraic(13)),
                   Conic_point_2(Algebraic(-3), Algebraic(9)));
  Conic_curve_2 c7(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(-3), Algebraic(9)),
                   Conic_point_2(Algebraic(0), Algebraic(0)));
  Conic_curve_2 c8(0,1,0,-1,0,0, CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(0), Algebraic(0)),
                   Conic_point_2(Algebraic(4), Algebraic(-2)));

  Conic_curve_2 c9(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                   Conic_point_2(Algebraic(-5), Algebraic(25)),
                   Conic_point_2(Algebraic(5), Algebraic(25)));
  Conic_curve_2 c10(58, 72, -48, 0, 0, -360);

  //This vector is used to store curves that will be used to create polycurve
  std::vector<Conic_curve_2> conic_curves;
  conic_curves.push_back(c9);

  //construct poly-curve
  Polycurve_conic_traits_2::Curve_2 conic_polycurve =
    construct_polycurve(conic_curves.begin(), conic_curves.end());

  Conic_curve_2 c11(0,1,0,-1,0,0,CGAL::COUNTERCLOCKWISE,
                     Conic_point_2(Algebraic(25), Algebraic(-5)),
                     Conic_point_2(Algebraic(0), Algebraic(0)));
  Conic_curve_2 c12(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
                     Conic_point_2(Algebraic(0), Algebraic(0)),
                     Conic_point_2(Algebraic(5), Algebraic(25)));
  conic_curves.clear();
  conic_curves.push_back(c11);
  conic_curves.push_back(c12);

  //construct poly-curve
  Polycurve_conic_traits_2::Curve_2 conic_polycurve_2 =
    construct_polycurve(conic_curves.begin(), conic_curves.end());

  /* VERY IMPORTANT
   * For efficiency reasons, we recommend users not to construct
   * x-monotone conic arc directly, but rather use the Make_x_monotone_2
   * functor supplied by the conic-arc traits class to convert conic curves
   * to x-monotone curves.
   */
  Conic_x_monotone_curve_2 xc3(c3);
  Conic_x_monotone_curve_2 xc4(c4);
  Conic_x_monotone_curve_2 xc5(c5);
  Conic_x_monotone_curve_2 xc6(c6);
  Conic_x_monotone_curve_2 xc7(c7);
  Conic_x_monotone_curve_2 xc8(c8);


  //This vector is used to store curves that will be used to create
  //X-monotone-polycurve
  std::vector<Conic_x_monotone_curve_2> xmono_conic_curves_2;
  xmono_conic_curves_2.push_back(xc5);
  xmono_conic_curves_2.push_back(xc3);
  xmono_conic_curves_2.push_back(xc4);


  //construct x-monotone poly-curve
  Pc_x_monotone_curve_2 conic_x_mono_polycurve_1 =
    construct_x_mono_polycurve(xmono_conic_curves_2.begin(),
                               xmono_conic_curves_2.end());

  xmono_conic_curves_2.clear();
  xmono_conic_curves_2.push_back(xc6);
  xmono_conic_curves_2.push_back(xc7);
  xmono_conic_curves_2.push_back(xc8);
  //construct x-monotone poly-curve
  Pc_x_monotone_curve_2 conic_x_mono_polycurve_2 =
    construct_x_mono_polycurve(xmono_conic_curves_2.begin(),
                               xmono_conic_curves_2.end());

  xmono_conic_curves_2.clear();
  xmono_conic_curves_2.push_back(xc5);

  Pc_x_monotone_curve_2 x_polycurve_push =
    construct_x_mono_polycurve(xmono_conic_curves_2.begin(),
                               xmono_conic_curves_2.end());
  Polycurve_conic_traits_2::X_monotone_subcurve_2 xcurve_push =
    Polycurve_conic_traits_2::X_monotone_subcurve_2(c5);
  //traits.construct_x_monotone_curve_2_object()(c5);

  xmono_conic_curves_2.clear();
  xmono_conic_curves_2.push_back(xc3);
  xmono_conic_curves_2.push_back(xc4);
  Pc_x_monotone_curve_2 base_curve =
    construct_x_mono_polycurve(xmono_conic_curves_2.begin(),
                               xmono_conic_curves_2.end());

  //curves for push_back
  Conic_curve_2 c13(1,1,0,-50,12,660,CGAL::COUNTERCLOCKWISE,
                    Conic_point_2(Algebraic(25), Algebraic(-7)),
                    Conic_point_2(Algebraic(25), Algebraic(-5)));
  Conic_curve_2 c14(0,1,0,-1,0,0,CGAL::COUNTERCLOCKWISE,
                    Conic_point_2(Algebraic(25), Algebraic(-5)),
                    Conic_point_2(Algebraic(0), Algebraic(0)));
  Conic_curve_2 c15(-1,0,0,0,1,0,CGAL::COUNTERCLOCKWISE,
                    Conic_point_2(Algebraic(0), Algebraic(0)),
                    Conic_point_2(Algebraic(5), Algebraic(25)));
  conic_curves.clear();
  conic_curves.push_back(c13);
  conic_curves.push_back(c14);
  Polycurve_conic_traits_2::Curve_2 base_curve_push_back =
    construct_polycurve(conic_curves.begin(), conic_curves.end());

  conic_curves.push_back(c15);
  Polycurve_conic_traits_2::Curve_2 Expected_push_back_result =
    construct_polycurve(conic_curves.begin(), conic_curves.end());

  // //checking the orientattion consistency
  // Conic_curve_2 c21(0,1,0,1,0,0,CGAL::CLOCKWISE,
  //                  Conic_point_2(Algebraic(9), Algebraic(-3)),
  //                  Conic_point_2(Algebraic(0), Algebraic(0)));
  // Conic_curve_2 c20(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
  //                   Conic_point_2(Algebraic(0), Algebraic(0)),
  //                   Conic_point_2(Algebraic(3), Algebraic(9)));
  //  Conic_x_monotone_curve_2 xc20(c20);
  //  Conic_x_monotone_curve_2 xc21(c21);
  //  xmono_conic_curves_2.clear();
  //  xmono_conic_curves_2.push_back(xc20);
  // xmono_conic_curves_2.push_back(xc21);
  // Pc_x_monotone_curve_2 eric_polycurve =
  //   construct_x_mono_polycurve(xmono_conic_curves_2.begin(),
  //                              xmono_conic_curves_2.end());
  // std::cout << "the polycurve is: " << eric_polycurve << std::endl;

  // std::cout<< std::endl;

  //check_compare_x_2(xc3, xc5);

  // check_equal();
  // std::cout<< std::endl;

   //check_intersect(conic_x_mono_polycurve_1, conic_x_mono_polycurve_2);
   //std::cout<< std::endl;

  // check_compare_end_points_xy_2();
  // std::cout<< std::endl;

  //check_split(conic_x_mono_polycurve_1, conic_x_mono_polycurve_2);
  // std::cout<< std::endl;

  //check_make_x_monotne_curve(conic_polycurve_2);
   //std::cout<< std::endl;

  // check_is_vertical();
  // std::cout<< std::endl;

  //check_compare_y_at_x_2();
  //std::cout<< std::endl;

  //adds the segment to the right.
  //check_push_back(base_curve_push_back, c15);
  //std::cout<< std::endl;

  //adds the segment to the left.
  //check_push_front(base_curve, xcurve_push);
  //std::cout<< std::endl;

  // check_are_mergable();
  // std::cout<< std::endl;

  // check_merge_2();
  // std::cout<< std::endl;

  // check_construct_opposite();
  // std::cout<< std::endl;

  // check_compare_y_at_x_right();
  // std::cout<< std::endl;

  // check_compare_y_at_x_left();
  // std::cout<< std::endl;
  //check_compare_points(conic_x_mono_polycurve_1);

  //number of segments
  //std::cout << "Number of segments: "
  //          << traits.number_of_points_2_object()(base_curve_push_back)
  //          << std::endl;

  check_trim(conic_x_mono_polycurve_1, atoi(argv[1]), atoi(argv[2]),
             atoi(argv[3]), atoi(argv[4]));
  std::cout << std::endl;

  //std::cout << (atoi(argv[1]) + atoi(argv[2])) << std::endl;
  // Conic_traits_2 con_traits;
  // Conic_curve_2 cc3(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE,
  //                   Conic_point_2(Algebraic(0), Algebraic(0)),
  //                   Conic_point_2(Algebraic(3), Algebraic(9)));
  // Conic_x_monotone_curve_2 xcc3(cc3);
  // Conic_point_2       ps2(0, 0);
  // Conic_point_2       pt2(3, 9);
  // std::cout << "conic curve is : " << xcc3 << std::endl;
  // Conic_x_monotone_curve_2 trimmed_curve =
  //   con_traits.trim_2_object()(xc3, ps2, pt2);
  // std::cout << "trimmed conic curve is : " << trimmed_curve << std::endl;

  return 0;
}

#endif
