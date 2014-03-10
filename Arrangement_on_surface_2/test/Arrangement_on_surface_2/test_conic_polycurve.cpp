// Testing the do_equal function

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
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>

//#include <CGAL/Arr_geometry_traits/Polyline_2.h>


////////////////////
//conic traits
////////////////////
typedef CGAL::CORE_algebraic_number_traits                                Nt_traits;
typedef Nt_traits::Rational                                               Rational;
typedef Nt_traits::Algebraic                                              Algebraic;
typedef CGAL::Cartesian<Rational>                                         Rat_kernel;
typedef Rat_kernel::Point_2                                               Rat_point_2;
typedef Rat_kernel::Segment_2                                             Rat_segment_2;
typedef Rat_kernel::Circle_2                                              Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                                        Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>       Conic_traits_2;
typedef Conic_traits_2::Point_2                                           Conic_point_2;
typedef Conic_traits_2::Curve_2                                           Conic_curve_2;
typedef Conic_traits_2::X_monotone_curve_2                                Conic_x_monotone_curve_2;
typedef CGAL::Arr_polyline_traits_2<Conic_traits_2>                       Polycurve_conic_traits_2;
//typedef Polycurve_conic_traits_2::Point_2                                 polypoint;


// typedef CGAL::Arr_polyline_traits_2<
//                                       CGAL::Arr_conic_traits_2<CGAL::Cartesian<BigRat>,
//                                       CGAL::Cartesian<Expr>, 
//                                       CGAL::CORE_algebraic_number_traits> 
//                                     >::Point_2   test_point_2;
  

 void check_equal()
 {
  
  bool are_equal;
  
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Equal_2 equal = traits.equal_2_object();
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();

  //create some curves
  Conic_point_2       ps1 (Rational(1,4), 4);
  Conic_point_2       pt1 (2, Rational(1,2));
  Conic_curve_2       c1 (0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE, ps1, pt1);
  

  Conic_point_2       ps2 (Rational(1,4), 4);
  Conic_point_2       pt2 (2, Rational(1,2));
  Conic_curve_2       c2 (0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE, ps2, pt2);


  Rat_point_2       ps3 (Rational(1,4), 4);
  Rat_point_2       pmid3(Rational(3,2), 2);
  Rat_point_2       pt3 (2, Rational(1,3));
  Conic_curve_2     c3 (ps3, pmid3, pt3);

  Rat_point_2       ps4 (1, 5);
  Rat_point_2       pmid4(Rational(3,2), 3);
  Rat_point_2       pt4 (3, Rational(1,3));
  Conic_curve_2     c4 (ps4, pmid4, pt4);


  // //make x_monotone
  Polycurve_conic_traits_2::X_monotone_curve_2 xmc1 = construct_x_monotone_curve_2(c1);
  Polycurve_conic_traits_2::X_monotone_curve_2 xmc2 = construct_x_monotone_curve_2(c2);
  Polycurve_conic_traits_2::X_monotone_curve_2 xmc3 = construct_x_monotone_curve_2(c3);
  Polycurve_conic_traits_2::X_monotone_curve_2 xmc4 = construct_x_monotone_curve_2(c4);  

  are_equal = equal(xmc1, xmc2);
  std::cout << "Two equal conic arcs are computed as:  " << ( (are_equal) ? "equal" : "Not equal") << std::endl;

  are_equal = equal(xmc3, xmc2);
  std::cout << "Two un-equal conic arcs are computed as:  " << ( (are_equal) ? "equal" : "Not equal") << std::endl;

  are_equal = equal(xmc3, xmc4);
  std::cout << "Two un-equal conic arcs are computed as:  " << ( (are_equal) ? "equal" : "Not equal") << std::endl;

 }

 void check_intersect()
 {

  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Intersect_2 intersect_2 = traits.intersect_2_object();
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();

  Rat_point_2       ps1 (Rational(1,4), 4);
  Rat_point_2       pmid1(Rational(3,2), 2);
  Rat_point_2       pt1 (2, Rational(1,3));
  Conic_curve_2     c1 (ps1, pmid1, pt1);

  Rat_point_2       ps2 (1, 10);
  Rat_point_2       pmid2(5, 4);
  Rat_point_2       pt2 (10, 1);
  Conic_curve_2     c2 (ps2, pmid2, pt2);

  Rat_point_2       ps3 (10, 1);
  Rat_point_2       pmid3(5, 4);
  Rat_point_2       pt3 (1, 10);
  Conic_curve_2     c3 (ps3, pmid3, pt3);


  //construct x-monotone curve (compatible with polyline class)
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 = construct_x_monotone_curve_2(c1);
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 = construct_x_monotone_curve_2(c2);
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc3 = construct_x_monotone_curve_2(c3);  


  std::vector<CGAL::Object> intersection_points;

  intersect_2(polyline_xmc1, polyline_xmc2, std::back_inserter(intersection_points));
  std::cout<< "For non-intersecting curves,  " << intersection_points.size() << " number of intersecting points were computed" << std::endl;

  intersection_points.clear();

  intersect_2(polyline_xmc2, polyline_xmc3, std::back_inserter(intersection_points));
  std::cout<< "For intersecting curves,  " << intersection_points.size() << " number of intersecting points were computed" << std::endl;

 }

void check_compare_end_points_xy_2()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Compare_endpoints_xy_2 compare_endpoints_xy_2 = traits.compare_endpoints_xy_2_object();

  //create some curves
  Conic_point_2       ps1 (Rational(1,4), 4);
  Conic_point_2       pt1 (2, Rational(1,2));
  Conic_curve_2       c1 (0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE, ps1, pt1);

  // Insert a parabolic arc that is supported by a parabola y = -x^2
  // (or: x^2 + y = 0) and whose endpoints are (-sqrt(3), -3) ~ (-1.73, -3)
  // and (sqrt(2), -2) ~ (1.41, -2). Notice that since the x-coordinates 
  // of the endpoints cannot be acccurately represented, we specify them
  // as the intersections of the parabola with the lines y = -3 and y = -2.
  // Note that the arc is clockwise oriented.
  Conic_curve_2   c2 = Conic_curve_2 (1, 0, 0, 0, 1, 0,       // The parabola.
                                      CGAL::CLOCKWISE,
                                       Conic_point_2 (-1.73, -3),    // Approximation of the source.
                                       0, 0, 0, 0, 1, 3,       // The line: y = -3.
                                       Conic_point_2 (1.41, -2),     // Approximation of the target.
                                       0, 0, 0, 0, 1, 2);      // The line: y = -2.
  CGAL_assertion (c2.is_valid());  

  //make polyline x-monotone curves
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 = construct_x_monotone_curve_2(c1);
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 = construct_x_monotone_curve_2(c2);

  CGAL::Comparison_result res = compare_endpoints_xy_2(polyline_xmc1);
  std::cout << "compare_end_points_xy_2 for counterclockwise curve: "<< (res == CGAL::SMALLER ? "SMALLER":
               (res == CGAL::LARGER ? "LARGER" : "EQUAL")) << std::endl;

  res = compare_endpoints_xy_2(polyline_xmc2);
  std::cout<< "compare_end_points_xy_2 for clockwise curve: "<< (res == CGAL::SMALLER ? "SMALLER":
               (res == CGAL::LARGER ? "LARGER" : "EQUAL")) << std::endl;
}

void check_split()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Split_2  split_2 = traits.split_2_object();

  //create a curve
  Rat_point_2       ps2 (1, 10);
  Rat_point_2       pmid2(5, 4);
  Rat_point_2       pt2 (10, 1);
  Conic_curve_2     c1 (ps2, pmid2, pt2);

  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 = construct_x_monotone_curve_2(c1);
  Polycurve_conic_traits_2::X_monotone_curve_2 split_curve_1, split_curve_2;
  Polycurve_conic_traits_2::Point_2 point_of_split = Polycurve_conic_traits_2::Point_2(5,4);

  split_2(polyline_xmc1, point_of_split, split_curve_1, split_curve_2);

  //TODO: check whether the split is correct  
}

void check_is_vertical()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Is_vertical_2 is_vertical = traits.is_vertical_2_object();

   //create a curve
  Rat_point_2       ps1 (1, 10);
  Rat_point_2       pmid1(5, 4);
  Rat_point_2       pt1 (10, 1);
  Conic_curve_2     c1 (ps1, pmid1, pt1);

  //make x-monotone curve
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 = construct_x_monotone_curve_2(c1);

  bool result = is_vertical(polyline_xmc1);
  std::cout << "Is_verticle:: Expected first result is not vertivle: Computed: " << ( (result)? "vertical" : "not vertical" ) << std::endl;

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
  
  //Basic_number_type x(rat_x), y(rat_y);
 // p = Conic_point_2( Rational(rat_x_num, rat_x_den), Rational(rat_y_num, rat_y_den) );
  //return true;
}

/*! */
template <typename stream>
bool read_orientation_and_end_points(stream& is, CGAL::Orientation& orient,
                                     Conic_point_2& source, Conic_point_2& target)
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

  std::cout << "line is: " << r << s << t << u << v << w << i_orient << r1 << s1 << t1 << u1 << v1 << w1 << r2 << s2 << t2 << u2 << v2 << w2 << std::endl;

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
  //std::cout<< "arc coefficients: " << r << " " << s << " " << t << " " << u << " " << v << " " << w << std::endl;
  //std::cout<< "Read Points : " << source.x() << " " << source.y() << " " << target.x() << " " << target.y() << std::endl;
  cv = Curve(r, s, t, u, v, w, orient, source, target);
  return true;
}
std::istream& skip_comments(std::istream& is,
                                                    std::string& line)
{
  while (std::getline(is, line))
    if (!line.empty() && (line[0] != '#')) break;
  return is;
}

bool check_compare_y_at_x_2()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Compare_y_at_x_2 cmp_y_at_x_2 = traits.compare_y_at_x_2_object();
  //polycurve constructors
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2 construct_x_mono_polycurve = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Construct_curve_2  construct_polycurve = traits.construct_curve_2_object();

   //create a curve
  Rat_point_2       ps1 (1, 10);
  Rat_point_2       pmid1(5, 4);
  Rat_point_2       pt1 (10, 1);
  Conic_curve_2     c1 (ps1, pmid1, pt1);

   //create a curve
  Rat_point_2       ps2 (10, 1);
  Rat_point_2       pmid2(15, 5);
  Rat_point_2       pt2 (20, 10);
  Conic_curve_2     c2 (ps2, pmid2, pt2);

  Conic_curve_2     c3(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE, Conic_point_2( Algebraic(0), Algebraic(0) ), Conic_point_2( Algebraic(3), Algebraic(9) ) );
  Conic_curve_2     c4(1,0,0,0,-1,0,CGAL::COUNTERCLOCKWISE, Conic_point_2( Algebraic(3), Algebraic(9) ), Conic_point_2( Algebraic(5), Algebraic(25) ) );

  std::vector<Conic_curve_2> conic_curves, conic_curves_2, Conic_curves_3;
  conic_curves.push_back(c1);
  conic_curves.push_back(c2);

  //conic_curves_2.push_back(c3);
  //conic_curves_2.push_back(c4);

  Conic_x_monotone_curve_2 xc1 (c1);
  Conic_x_monotone_curve_2 xc2 (c2);
  Conic_x_monotone_curve_2 xc3 (c3);
  Conic_x_monotone_curve_2 xc4 (c4);

  
  std::vector<Conic_x_monotone_curve_2> xmono_conic_curves, xmono_conic_curves_2;
  /* VERY IMPORTANT
  For efficiency reasons, we recommend users not to construct x-monotone conic arc directly, 
  but rather use the Make_x_monotone_2 functor supplied by the conic-arc traits class to convert conic curves to x-monotone curves.
  */
  xmono_conic_curves.push_back(xc1);
  xmono_conic_curves.push_back(xc2);
  xmono_conic_curves_2.push_back(xc3);
  xmono_conic_curves_2.push_back(xc4);

  ////////////////////
  //Reading from a file
  ////////////////
  // std::string filename = "data/polycurves_conics/compare_y_at_x.xcv";
  // std::ifstream cv_stream(filename);
  
  // if (!cv_stream.is_open()) 
  // {
  //   std::cerr << "Cannot open file " << filename << "!" << std::endl;
  //   return false;
  // }

  // std::string line;
  
  // while (skip_comments(cv_stream, line)) 
  // {
  //   std::cout<< "line fiele is :" << line << std::endl;
  //   std::istringstream line_stream(line);
  //   Conic_curve_2 cv;
  //   read_general_arc(line_stream, cv);
  //   //Conic_curves_3.push_back(cv);
  //   line_stream.clear();
  // }
  
  // cv_stream.close();

  ////////
  ///////




  //construct x-monotone poly-curve
  Polycurve_conic_traits_2::X_monotone_curve_2 conic_x_mono_polycurve = construct_x_mono_polycurve(xmono_conic_curves.begin(), xmono_conic_curves.end());
  Polycurve_conic_traits_2::X_monotone_curve_2 conic_x_mono_polycurve_2 = construct_x_mono_polycurve(xmono_conic_curves_2.begin(), xmono_conic_curves_2.end());

  //construct poly-curve
  Polycurve_conic_traits_2::Curve_2 conic_polycurve = construct_polycurve( conic_curves.begin(), conic_curves.end() );
  //Polycurve_conic_traits_2::Curve_2 conic_polycurve_2 = construct_polycurve( conic_curves_2.begin(), conic_curves_2.end() );

  //make x-monotone curve
  //Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 = construct_x_monotone_curve_2(c1);
 
  //create points
  Polycurve_conic_traits_2::Point_2 point_above_line = Polycurve_conic_traits_2::Point_2(2,10), 
                                    point_below_line = Polycurve_conic_traits_2::Point_2(4,7), 
                                    point_on_line    = Polycurve_conic_traits_2::Point_2(2,4);

  CGAL::Comparison_result result;

  result =  cmp_y_at_x_2(point_above_line, conic_x_mono_polycurve_2);
  std::cout << "Compare_y_at_x_2:: for point above the curve computed Answer is:  "<< (result == CGAL::SMALLER ? "Below":
               (result == CGAL::LARGER ? "Above" : "On-line")) << std::endl;

  result =  cmp_y_at_x_2(point_below_line, conic_x_mono_polycurve_2);
  std::cout << "Compare_y_at_x_2:: for point below the curve computed Answer is:  "<< (result == CGAL::SMALLER ? "Below":
               (result == CGAL::LARGER ? "Above" : "On-line")) << std::endl;

  result =  cmp_y_at_x_2(point_on_line, conic_x_mono_polycurve_2);
  std::cout << "Compare_y_at_x_2:: for point on the curve computed Answer is:  "<< (result == CGAL::SMALLER ? "Below":
               (result == CGAL::LARGER ? "Above" : "On-line")) << std::endl;

  return true;
}

void check_push_back()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Push_back_2  push_back_2 = traits.push_back_2_object();

  Polycurve_conic_traits_2::Curve_2 polycurve;
  Polycurve_conic_traits_2::X_monotone_curve_2 x_monotone_polycurve;

   //create a curve
  Rat_point_2       ps1 (1, 10);
  Rat_point_2       pmid1(5, 4);
  Rat_point_2       pt1 (10, 1);
  Conic_curve_2     c1 (ps1, pmid1, pt1);

  Rat_point_2       ps2 (10, 1);
  Rat_point_2       pmid2(15, 14);
  Rat_point_2       pt2 (20, 10);
  Conic_curve_2     c2 (ps2, pmid2, pt2);
  

  // push_back_2(polycurve, c1);
  // std::cout<< "Push_back_2::size of polycurve after 1 segment push_back: " << polycurve.size() << std::endl;

  // push_back_2(polycurve, c2);
  // std::cout<< "Push_back_2::size of polycurve after 2 segment push_back: " << polycurve.size() << std::endl;

  // push_back_2( x_monotone_polycurve, Polycurve_conic_traits_2::X_monotone_segment_2(c1) );
  // std::cout<< "Push_back_2::size of x-monotone polycurve after 1 x-monotone segment push_back: " << x_monotone_polycurve.size() << std::endl;
  
  //error
  //push_back_2( x_monotone_polycurve, Polycurve_conic_traits_2::X_monotone_segment_2(c2) );
  //std::cout<< "Push_back_2::size of x-monotone polycurve after 2 x-monotone segment push_back: " << x_monotone_polycurve.size() << std::endl;

}

void check_push_front()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Push_front_2  push_front_2 = traits.push_front_2_object();

  Polycurve_conic_traits_2::Curve_2 polycurve;
  Polycurve_conic_traits_2::X_monotone_curve_2 x_monotone_polycurve;

   //create a curve
  Rat_point_2       ps1 (1, 10);
  Rat_point_2       pmid1(5, 4);
  Rat_point_2       pt1 (10, 1);
  Conic_curve_2     c1 (ps1, pmid1, pt1);

  Rat_point_2       ps2 (10, 1);
  Rat_point_2       pmid2(15, 14);
  Rat_point_2       pt2 (20, 10);
  Conic_curve_2     c2 (ps2, pmid2, pt2);

  // push_front_2(polycurve, c1);
  // std::cout<< "Push_front_2::size of polycurve after 1 segment push_back: " << polycurve.size() << std::endl;

  // push_front_2(polycurve, c2);
  // std::cout<< "Push_front_2::size of polycurve after 2 segment push_back: " << polycurve.size() << std::endl;

  // push_front_2( x_monotone_polycurve, Polycurve_conic_traits_2::X_monotone_segment_2(c1) );
  // std::cout<< "Push_front_2::size of x-monotone polycurve after 1 x-monotone segment push_back: " << x_monotone_polycurve.size() << std::endl;
  
  //error
  //push_front_2( x_monotone_polycurve, Polycurve_conic_traits_2::X_monotone_segment_2(c2) );
  //std::cout<< "Push_front_2::size of x-monotone polycurve after 2 x-monotone segment push_back: " << x_monotone_polycurve.size() << std::endl;
}

void check_are_mergable()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Are_mergeable_2  are_mergeable_2 = traits.are_mergeable_2_object();

   //create a curve
  Rat_point_2       ps1 (1, 10);
  Rat_point_2       pmid1(5, 4);
  Rat_point_2       pt1 (10, 1);
  Conic_curve_2     c1 (ps1, pmid1, pt1);

  Rat_point_2       ps2 (10, 1);
  Rat_point_2       pmid2(15, 14);
  Rat_point_2       pt2 (20, 20);
  Conic_curve_2     c2 (ps2, pmid2, pt2);

  Rat_point_2       ps3 (Rational(1,4), 4);
  Rat_point_2       pmid3(Rational(3,2), 2);
  Rat_point_2       pt3 (2, Rational(1,3));
  Conic_curve_2     c3 (ps3, pmid3, pt3);

  //construct x-monotone curve (compatible with polyline class)
   Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 = construct_x_monotone_curve_2(c1);
   Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 = construct_x_monotone_curve_2(c2);
   Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc3 = construct_x_monotone_curve_2(c3); 

  bool result;

  result = are_mergeable_2(polyline_xmc1, polyline_xmc2);
  std::cout << "Are_mergable:: Mergable x-monotone polycurves are Computed as: " << ( (result)? "Mergable" : "Not-Mergable" ) << std::endl;

  result = are_mergeable_2(polyline_xmc1, polyline_xmc3);
  std::cout << "Are_mergable:: Non-Mergable x-monotone polycurves are Computed as: " << ( (result)? "Mergable" : "Not-Mergable" ) << std::endl;

}

void check_merge_2()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Merge_2  merge_2 = traits.merge_2_object();

  //create a curve
  Rat_point_2       ps1 (1, 10);
  Rat_point_2       pmid1(5, 4);
  Rat_point_2       pt1 (10, 1);
  Conic_curve_2     c1 (ps1, pmid1, pt1);

  Rat_point_2       ps2 (10, 1);
  Rat_point_2       pmid2(15, 14);
  Rat_point_2       pt2 (20, 20);
  Conic_curve_2     c2 (ps2, pmid2, pt2);

//construct x-monotone curve (compatible with polyline class)
 Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 = construct_x_monotone_curve_2(c1);
 Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 = construct_x_monotone_curve_2(c2);
 
 Polycurve_conic_traits_2::X_monotone_curve_2 merged_xmc; 


 merge_2(polyline_xmc1, polyline_xmc2, merged_xmc);
 std::cout<< "Merge_2:: Mergable x-monotone curves merged successfully"<< std:: endl;

}

void check_construct_opposite()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Construct_opposite_2  construct_opposite_2 = traits.construct_opposite_2_object();

  //create a curve
  Rat_point_2       ps1 (1, 10);
  Rat_point_2       pmid1(5, 4);
  Rat_point_2       pt1 (10, 1);
  Conic_curve_2     c1 (ps1, pmid1, pt1);

  //construct x-monotone curve (compatible with polyline class)
 Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 = construct_x_monotone_curve_2(c1);
 Polycurve_conic_traits_2::X_monotone_curve_2 polyline_opposite_curve = construct_opposite_2(polyline_xmc1);

 std::cout<< "Construct_opposite_2:: Opposite curve created";


}

void check_compare_y_at_x_right()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Compare_y_at_x_right_2 cmp_y_at_x_right_2 = traits.compare_y_at_x_right_2_object();

  //create constructing curves
  Rat_point_2       ps2 (1, 10);
  Rat_point_2       pmid2(5, 4);
  Rat_point_2       pt2 (10, 1);
  Conic_curve_2     c1 (ps2, pmid2, pt2);

  Rat_point_2       ps3 (10, 1);
  Rat_point_2       pmid3(5, 4);
  Rat_point_2       pt3 (1, 10);
  Conic_curve_2     c2 (ps3, pmid3, pt3);


  //construct x-monotone curve (compatible with polyline class)
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 = construct_x_monotone_curve_2(c1);
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 = construct_x_monotone_curve_2(c2); 
  Polycurve_conic_traits_2::Point_2 intersection_point = Polycurve_conic_traits_2::Point_2(5,4);

  CGAL::Comparison_result result;

  result = cmp_y_at_x_right_2(polyline_xmc1, polyline_xmc2, intersection_point);
  std::cout << "Compare_y_at_x_right:: Expected Answer: equal, Computed answer:  "<< (result == CGAL::SMALLER ? "smaller":
               (result == CGAL::LARGER ? "Larger" : "equal")) << std::endl;
}

void check_compare_y_at_x_left()
{
  Polycurve_conic_traits_2 traits;
  Polycurve_conic_traits_2::Construct_x_monotone_curve_2  construct_x_monotone_curve_2 = traits.construct_x_monotone_curve_2_object();
  Polycurve_conic_traits_2::Compare_y_at_x_left_2 cmp_y_at_x_left_2 = traits.compare_y_at_x_left_2_object();

  //create constructing curves
  Rat_point_2       ps2 (1, 10);
  Rat_point_2       pmid2(5, 4);
  Rat_point_2       pt2 (10, 1);
  Conic_curve_2     c1 (ps2, pmid2, pt2);

  Rat_point_2       ps3 (10, 1);
  Rat_point_2       pmid3(5, 4);
  Rat_point_2       pt3 (1, 10);
  Conic_curve_2     c2 (ps3, pmid3, pt3);


  //construct x-monotone curve (compatible with polyline class)
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc1 = construct_x_monotone_curve_2(c1);
  Polycurve_conic_traits_2::X_monotone_curve_2 polyline_xmc2 = construct_x_monotone_curve_2(c2); 
  Polycurve_conic_traits_2::Point_2 intersection_point = Polycurve_conic_traits_2::Point_2(5,4);

  CGAL::Comparison_result result;

  result = cmp_y_at_x_left_2(polyline_xmc1, polyline_xmc2, intersection_point);
  std::cout << "Compare_y_at_x_left:: Expected Answer: equal, Computed answer:  "<< (result == CGAL::SMALLER ? "smaller":
               (result == CGAL::LARGER ? "Larger" : "equal")) << std::endl;
}

int main ()
{

  // std::cout<< std::endl;

  // check_equal();
  // std::cout<< std::endl;
  
  // check_intersect();
  // std::cout<< std::endl;
  
  // check_compare_end_points_xy_2();
  // std::cout<< std::endl;

  // check_split();
  // std::cout<< std::endl;

  // check_is_vertical();
  // std::cout<< std::endl;

   check_compare_y_at_x_2();
   std::cout<< std::endl;

  // check_push_back();
  // std::cout<< std::endl;

  // check_push_front();
  // std::cout<< std::endl;

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

  return 0;
}

#endif