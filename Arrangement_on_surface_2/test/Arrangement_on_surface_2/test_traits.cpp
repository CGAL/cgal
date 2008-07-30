#include <CGAL/basic.h>
#include "test_configuration.h"
#include <iostream>

#if (TEST_TRAITS == CORE_CONIC_TRAITS) && !defined(CGAL_USE_CORE)

int main ()
{
  bool   UNTESTED_TRAITS_AS_CORE_IS_NOT_ISTALLED;
  std::cout << std::endl
            << "WARNING: Core is not installed, "
            << "skipping the test of the conic traits ..."
            << std::endl; 
  return (0);
}
#else

#include <CGAL/assertions.h>
#include <CGAL/Arrangement_2.h>

#include <vector>

#include "test_traits.h"
#include "Traits_test.h"

// Arrangement types:
typedef CGAL::Arr_default_dcel<Traits>                  Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel>               Arr;

// Traits types:
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef Traits::Curve_2                                 Curve_2;

int main (int argc, char * argv[])
{
  CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
  CGAL::set_warning_behaviour(CGAL::THROW_EXCEPTION);
  Traits_test<Traits> test(argc, argv);
  bool rc  = test.start();
  return (rc) ? 0 : -1;
}

#if TEST_TRAITS == CORE_CONIC_TRAITS || \
    TEST_TRAITS == RATIONAL_ARC_TRAITS

// conic traits and rational traits use same number 
// type CORE:Expr so this code can be shared

/*! Read a point */

template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_point(stream & is, Point_2 & p)
{
  Rational rat_x,rat_y;
  is >> rat_x >> rat_y;
  Basic_number_type x(rat_x), y(rat_y);
  p = Point_2(x, y);
  return true;
}

#endif

#if TEST_TRAITS == SPHERICAL_ARC_TRAITS

/*! Read a point */

template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_point(stream & is, Point_2 & p)
{
  Basic_number_type x, y, z;
  is >> x >> y >> z;
  p = Point_2(x, y, z);
  return true;
}

/*! Read a xcurve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream & is, X_monotone_curve_2 & xcv)
{
  Point_2 p1,p2;
  read_point(is, p1);
  read_point(is, p2);
  CGAL_assertion(p1 != p2);
  xcv = X_monotone_curve_2(p1, p2);
  return true;
}

/*! Read a curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_curve(stream & is, Curve_2 & cv)
{
  Point_2 p1, p2;
  read_point(is, p1);
  read_point(is, p2);
  CGAL_assertion(p1 != p2);
  cv = Curve_2(p1, p2);
  return true;
}

#endif

#if TEST_TRAITS == RATIONAL_ARC_TRAITS

/*! Read a xcurve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream & is, X_monotone_curve_2 & xcv)
{
  Curve_2 tmp_cv;
  if (!read_curve(is,tmp_cv))
    return false;
  xcv=X_monotone_curve_2(tmp_cv);
  return true;
}

template <class stream>
bool read_coefficients(stream & is,Rat_vector & coeffs)
{
  unsigned int num_coeffs;
  Rational rat;
  is >> num_coeffs;
  coeffs.clear();
  for (unsigned int j = 0; j < num_coeffs; j++) {
    is >> rat;
    coeffs.push_back(rat);
  }
  return true;
}

/*! Read a curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::read_curve(stream & is, Curve_2 & cv)
{
  // Get the arc type:
  Rat_vector p_coeffs, q_coeffs;
  Algebraic src, trg;
  int dir;
  char type;
  is >> type;
  if (type == 'a' || type == 'A') 
  {
    //Default constructor
    cv = Curve_2();
    return true;
  }
  else if (type == 'b' || type == 'B') 
  {
    //Constructor of a whole polynomial curve
    if (read_coefficients(is,p_coeffs))
      cv = Curve_2(p_coeffs);
    else
      return false;
    return true;
  }
  else if (type == 'c' || type == 'C') 
  {
    //Constructor of a polynomial ray
    if (!read_coefficients(is,p_coeffs))
      return false;
    is >> src >> dir;
    cv = Curve_2(p_coeffs ,src ,(dir == 0 ? false : true));
    return true;
  }
  else if (type == 'd' || type == 'D') 
  {
    //Constructor of a polynomial arc
    if (!read_coefficients(is,p_coeffs))
      return false;
    is >> src >> trg;
    cv = Curve_2(p_coeffs ,src ,trg);
    return true;
  }
  else if (type == 'e' || type == 'E') 
  {
    //Constructor of a whole rational function
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_coefficients(is,q_coeffs))
      return false;
    cv = Curve_2(p_coeffs ,q_coeffs);
    return true;
  }
  else if (type == 'f' || type == 'F') 
  {
    //Constructor of a ray of a rational function
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_coefficients(is,q_coeffs))
      return false;
    is >> src >> dir;
    cv = Curve_2(p_coeffs ,q_coeffs ,src ,(dir == 0 ? false : true));
    return true;
  }
  else if (type == 'g' || type == 'G') 
  {
    //Constructor of a bounded rational arc
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_coefficients(is,q_coeffs))
      return false;
    is >> src >> trg;
    cv = Curve_2(p_coeffs ,q_coeffs ,src ,trg);
    return true;
  }
  // If we reached here, we have an unknown rational arc type:
  std::cerr << "Illegal rational arc type specification: " << type << "."
            << std::endl;
  return (false);
}

#endif

#if TEST_TRAITS == POLYLINE_TRAITS || TEST_TRAITS == NON_CACHING_POLYLINE_TRAITS

template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_xcurve(stream & is,
            Traits::X_monotone_curve_2 & xcv)
{
  unsigned int num_points;
  is >> num_points;
  std::vector<Point_2> points;
  points.clear();
  for (unsigned int j = 0; j < num_points; j++) {
    Basic_number_type x, y;
    is >> x >> y;
    Point_2 p(x, y);
    points.push_back(p);
  }
  xcv =
    CGAL::
    Arr_polyline_traits_2<Segment_traits>::X_monotone_curve_2(points.begin(),
                                                              points.end());
  return true;
}

template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_curve(stream & is,
           Traits::Curve_2 & cv)
{
  unsigned int num_points;
  is >> num_points;
  std::vector<Point_2> points;
  points.clear();
  for (unsigned int j = 0; j < num_points; j++) {
    Basic_number_type x, y;
    is >> x >> y;
    Point_2 p(x, y);
    points.push_back(p);
  }
  cv =
    CGAL::
    Arr_polyline_traits_2<Segment_traits>::Curve_2(points.begin(),
                                                   points.end());
  return true;
}

#elif TEST_TRAITS == LINEAR_TRAITS

template <>
template <class stream>
bool
Traits_base_test<Traits >::read_xcurve(stream & is, X_monotone_curve_2 & xcv)
{
  is >> xcv;
  return true;
}

template <>
template <class stream>
bool
Traits_base_test<Traits >::read_curve(stream & is, Curve_2 & cv)
{
  is >> cv;
  return true;
}

#elif TEST_TRAITS == CORE_CONIC_TRAITS

/*! */
template <class stream>
bool read_orientation(stream & is, CGAL::Orientation & orient)
{
  int i_orient;
  is >> i_orient;
  orient = (i_orient > 0) ? CGAL::COUNTERCLOCKWISE :
    (i_orient < 0) ? CGAL::CLOCKWISE : CGAL::COLLINEAR;
  return true;
}

/*! */
template <class stream>
bool read_app_point(stream & is, Point_2 & p)
{
  double x, y;
  is >> x >> y;
  p = Point_2(Algebraic(x), Algebraic(y));
  return true;
}

/*! */
template <class stream>
bool read_orientation_and_end_points(stream & is, CGAL::Orientation & orient,
                                     Point_2 & source, Point_2 & target)
{
  // Read the orientation.
  if (!read_orientation(is, orient)) return false;

  // Read the end points of the arc and create it.
  if (!read_app_point(is, source)) return false;
  if (!read_app_point(is, target)) return false;
  return true;
}

/*! Read a circle or an ellipse */
template <class stream>
bool read_ellipse(stream & is, bool & is_circle, Rat_circle & circle,
                  Rational & r, Rational & s,
                  Rational & t, Rational & u, Rational & v, Rational & w)
{
  // Read the ellipse (using the format "a b x0 y0"):
  //              2               2
  //   ( x - x0 )      ( y - y0 )
  //    --------   +   --------   = 1
  //       a              b
  //
  Rational a, b, x0, y0;
  is >> a >> b >> x0 >> y0;

  if (a == b) {
    is_circle = true;
    circle = Rat_circle (Rat_point (x0, y0), a);
  }
  else {
    r = 1/a;
    s = 1/b;
    t = 0;
    u = -2*x0*b;
    v = -2*y0*a;
    w = x0*x0*b + y0*y0*a - a*b;
  }
  return true;
}

/*! */
template <class stream, class Curve>
bool read_partial_ellipse(stream & is, Curve & cv)
{
  bool is_circle;               // Is this a circle.
  Rat_circle circle;
  Rational r, s, t, u, v, w;
  if (!read_ellipse(is, is_circle, circle, r, s, t, u, v, w)) return false;
  CGAL::Orientation orient;
  Point_2 source, target;
  if (!read_orientation_and_end_points(is, orient, source, target))  //!!!
    return false;

  // Create the conic (or circular) arc.
  cv = (is_circle) ? Curve (circle, orient, source, target) :
    Curve (r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! */
template <class stream, class Curve>
bool read_full_ellipse(stream & is, Curve & cv)
{
  bool is_circle;               // Is this a circle.
  Rat_circle circle;
  Rational r, s, t, u, v, w;
  if (!read_ellipse(is, is_circle, circle, r, s, t, u, v, w))
    return false;

  // Create a full ellipse (or circle).
  cv = (is_circle) ? Curve (circle) : Curve (r, s, t, u, v, w);
  return true;
}

/*! Read a hyperbola */
template <class stream, class Curve>
bool read_hyperbola(stream & is, Curve & cv)
{
  // Read the hyperbola (using the format "a b x0 y0"):
  //              2              2
  //    ( x - x0 )     ( y - y0 )
  //     --------   -   --------   = 1
  //       a               b
  //
  Rational a, b, x0, y0;
  is >> a >> b >> x0 >> y0;

  Rational r = b;
  Rational s= -a;
  Rational t = 0;
  Rational u = -2*x0*b;
  Rational v = 2*y0*a;
  Rational w = x0*x0*b - y0*y0*a - a*b;

  CGAL::Orientation orient;
  Point_2 source, target;
  if (!read_orientation_and_end_points(is, orient, source, target)) //!!!
    return false;

  // Create the conic (or circular) arc.
  cv = Curve (r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! Read a hyperbola */
template <class stream, class Curve>
bool read_parabola(stream & is, Curve & cv)
{
  // Read the parabola (using the format "c x0 y0"):
  //
  //          2
  //  (x - x0) = 4c*(y - y0)
  //
  Rational c, x0, y0;
  is >> c >> x0 >> y0;

  Rational r = 1;
  Rational s = 0;
  Rational t = 0;
  Rational u = -2*x0;
  Rational v = -4*c;
  Rational w = x0*x0 + 4*c*y0;

  CGAL::Orientation orient;
  Point_2 source, target;
  if (!read_orientation_and_end_points(is, orient, source, target))
    return false;

  // Create the conic (or circular) arc.
  cv = Curve (r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! */
template <class stream, class Curve>
bool read_segment(stream & is, Curve & cv)
{

  // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
  Rational x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  
  Rat_point source (x1, y1);
  Rat_point target (x2, y2);
  Rat_segment segment (source, target);

  // Create the segment.
  cv = Curve (segment);
  return true;
}

/*! */
template <class stream, class Curve>
bool read_general_arc(stream & is, Curve & cv)
{
  // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
  Rational r, s, t, u, v, w;                // The conic coefficients.
  is >> r >> s >> t >> u >> v >> w;
    // Read the orientation.
  int i_orient;
  is >> i_orient;
  CGAL::Orientation orient = (i_orient > 0) ? CGAL::COUNTERCLOCKWISE :
    (i_orient < 0) ? CGAL::CLOCKWISE : CGAL::COLLINEAR;

  // Read the approximated source, along with a general conic 
  // <r_1,s_1,t_1,u_1,v_1,w_1> whose intersection with <r,s,t,u,v,w>
  // defines the source.
  Point_2 app_source;
  if (!read_app_point(is, app_source)) return false;
  Rational r1, s1, t1, u1, v1, w1;
  is >> r1 >> s1 >> t1 >> u1 >> v1 >> w1;

  // Read the approximated target, along with a general conic 
  // <r_2,s_2,t_2,u_2,v_2,w_2> whose intersection with <r,s,t,u,v,w>
  // defines the target.
  Point_2 app_target;
  if (!read_app_point(is, app_target)) return false;

  Rational r2, s2, t2, u2, v2, w2;
  is >> r2 >> s2 >> t2 >> u2 >> v2 >> w2;

  // Create the conic arc.
  cv = Curve (r, s, t, u, v ,w, orient,
              app_source, r1, s1, t1, u1, v1, w1,
              app_target, r2, s2, t2, u2, v2, w2);
  return true;
}

/*! */
template <class stream, class Curve>
bool read_general_curve(stream & is, Curve & cv)
{
  Rational r, s, t, u, v, w;                // The conic coefficients.
  // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
  is >> r >> s >> t >> u >> v >> w;
  CGAL::Orientation orient;
  Point_2 source, target;
  if (!read_orientation_and_end_points(is, orient, source, target))
    return false;

  // Create the conic (or circular) arc.
  cv = Curve (r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! Read an x-monotone conic curve */
template <>
template <class stream>
bool
Traits_base_test<Traits>::
read_xcurve(stream & is,
            X_monotone_curve_2 & xcv)
{
  Curve_2 tmp_cv;
  if (!read_curve(is,tmp_cv))
    return false;
  xcv=X_monotone_curve_2(tmp_cv);
  return true;
}

/*! Read a general conic curve */
template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_curve(stream & is,
           Curve_2 & cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'f' || type == 'F') 
  {
    return read_full_ellipse(is, cv);
  }
  else if (type == 's' || type == 'S') 
  {
    return read_segment(is, cv);
  }
  else if (type == 'i' || type == 'I') 
  {
    return read_general_arc(is, cv);
  }
  else if (type == 'c' || type == 'C') 
  {
    // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
    Rational r, s, t, u, v, w;
    is >> r >> s >> t >> u >> v >> w;
    // Create a full conic (should work only for ellipses).
    cv = Curve_2(r, s, t, u, v, w);
    return (true);
  }
  else if (type == 'e' || type == 'E') 
  {
    return read_partial_ellipse(is, cv);
  }
  else if (type == 'h' || type == 'H') 
  {
    return read_hyperbola(is, cv);
  }
  else if (type == 'p' || type == 'P') 
  {
    return read_parabola(is, cv);
  }
  else if (type == 'a' || type == 'A') 
  {
    return read_general_curve(is, cv);
  }

  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal conic type specification: " << type << "."
	    << std::endl;
  return (false);
}

#elif TEST_TRAITS == CIRCLE_SEGMENT_TRAITS

template <class stream>
bool read_ort_point(stream & is, Point_2 & p)
{
  bool is_rat;
  typename Point_2::CoordNT ort_x , ort_y;
  Number_type alpha,beta,gamma;
  is >> is_rat;
  if (is_rat)
  {
    is >> alpha;
    ort_x=Point_2::CoordNT(alpha);
  }
  else
  {
    is >> alpha >> beta >> gamma;
    ort_x=Point_2::CoordNT(alpha,beta,gamma);
  }
  is >> is_rat;
  if (is_rat)
  {
    is >> alpha;
    ort_y=Point_2::CoordNT(alpha);
  }
  else
  {
    is >> alpha >> beta >> gamma;
    ort_y=Point_2::CoordNT(alpha,beta,gamma);
  }
  p = Point_2(ort_x , ort_y);
  return true;
}

/*! Read an x-monotone circle segment curve */
template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_xcurve(stream & is,X_monotone_curve_2 & xcv)
{
  bool ans=true;
  char type;
  is >> type;
  if (type == 'z' || type == 'Z')
  {
    Line_2 l;
    Point_2 ps,pt;
    is >> l;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    xcv=X_monotone_curve_2(l,ps,pt);
    return ans;
  }
  else if (type == 'y' || type == 'Y')
  {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    xcv=X_monotone_curve_2(ps,pt);
    return true;
  }
  else if (type == 'x' || type == 'X')
  {
    Circle_2 c;
    Point_2 ps,pt;
    is >> c;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    xcv=X_monotone_curve_2(c,ps,pt,c.orientation());
    return ans;
  }
  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal circle segment type specification: " << type << std::endl;
  return false;
}

/*! Read a general circle segment curve */
template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_curve(stream & is,Curve_2 & cv)
{
  bool ans=true;
  char type;
  is >> type;
  if (type == 'a' || type == 'A')
  {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    Segment_2 s(ps,pt);
    cv=Curve_2(s);
    return true;
  }
  else if (type == 'b' || type == 'B')
  {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    cv=Curve_2(ps,pt);
    return true;
  }
  else if (type == 'c' || type == 'C')
  {
    Line_2 l;
    Point_2 ps,pt;
    is >> l;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv=Curve_2(l,ps,pt);
    return ans;
  }
  else if (type == 'd' || type == 'D')
  {
    Circle_2 c;
    is >> c;
    cv=Curve_2(c);
    return true;
  }
  else if (type == 'e' || type == 'E')
  {
    Rat_point_2 p;
    Rat_nt r;
    int orient;
    is >> p >> r >> orient;
    cv=Curve_2(p,r,static_cast<CGAL::Orientation>(orient));
    return true;
  }
  else if (type == 'f' || type == 'F')
  {
    Circle_2 c;
    Point_2 ps,pt;
    is >> c;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv=Curve_2(c,ps,pt);
    return ans;
  }
  else if (type == 'g' || type == 'G')
  {
    Rat_point_2 p;
    Rat_nt r;
    int orient;
    Point_2 ps,pt;
    is >> p >> r >> orient;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv=Curve_2(p,r,static_cast<CGAL::Orientation>(orient),ps,pt);
    return ans;
  }
  else if (type == 'h' || type == 'H')
  {
    Rat_point_2 ps,pm,pt;
    is >> ps >> pm >> pt;
    cv=Curve_2(ps,pm,pt);
    return true;
  }
  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal circle segment type specification: " << type << std::endl;
  return false;
}

#elif TEST_TRAITS == BEZIER_TRAITS

template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_point(stream & is, Point_2 & p)
{
  Rational rat_x,rat_y;
  is >> rat_x >> rat_y;
  p = Point_2(rat_x, rat_y);
  return true;
}

/*! Read an x-monotone bezier curve */

template <>
template <class stream>
bool
Traits_base_test<Traits>::read_xcurve(stream & is, X_monotone_curve_2 & xcv)
{
  std::list<CGAL::Object>                  x_objs;
  std::list<CGAL::Object>::const_iterator  xoit;
  Curve_2 tmp_cv;
  is >> tmp_cv;
  Rational B_psx = Rational(tmp_cv.control_point(0).x());
  Rational B_psy = Rational(tmp_cv.control_point(0).y());
  Rational B_ptx = Rational(tmp_cv.control_point(tmp_cv.number_of_control_points()-1).x());
  Rational B_pty = Rational(tmp_cv.control_point(tmp_cv.number_of_control_points()-1).y());
  Point_2 B_ps(B_psx, B_psy);
  Point_2 B_pt(B_ptx, B_pty);
  Traits::Make_x_monotone_2 make_x_monotone = this->m_traits.make_x_monotone_2_object();
  make_x_monotone (tmp_cv, std::front_inserter (x_objs));
  xoit = x_objs.begin();
  if (CGAL::assign (xcv, *xoit))
    return true;
  return false;
}

/*! Read a general bezier curve */
template <>
template <class stream>
bool
Traits_base_test<Traits >::read_curve(stream & is,Curve_2 & cv)
{
  is >> cv;
  return true;
}

#elif TEST_TRAITS == LINE_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

/*! Read an arc point */
template < typename T_Traits , typename stream >
bool read_arc_point(stream & is, Point_2 & p)
{
  Basic_number_type x, y;
  is >> x >> y;
  Circular_kernel::Point_2 lp(x, y);
  p = Point_2(lp);
  return true;
}

bool is_deg_1(char c)
{
  return (c=='z' || c=='Z') || (c=='y' || c=='Y') || (c=='x' || c=='X') ||
         (c=='w' || c=='W') || (c=='v' || c=='V') || (c=='l' || c=='L');
}

bool is_deg_2(char c)
{
  return (c=='b' || c=='B') || (c=='c' || c=='C') ||
         (c=='d' || c=='D') || (c=='e' || c=='E');
}

#if TEST_TRAITS == LINE_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS
template <class stream>
Circular_kernel::Line_arc_2 read_line(char type,stream & is)
{
  if (type == 'z' || type == 'Z')
  {
    Circular_kernel::Line_2 l_temp;
    Circular_kernel::Circle_2 c_temp1,c_temp2;
    bool b1,b2;
    is >> l_temp >> c_temp1 >> b1 >> c_temp2 >> b2;
    return Circular_kernel::Line_arc_2(l_temp,c_temp1,b1,c_temp2,b2);
  }
  else if (type == 'y' || type == 'Y')
  {
    Circular_kernel::Line_2 l_temp,l_temp1,l_temp2;
    is >> l_temp >> l_temp1 >> l_temp2;
    return Circular_kernel::Line_arc_2(l_temp,l_temp1,l_temp2);
  }
  else if (type == 'x' || type == 'X')
  {
    Circular_kernel::Line_2 l_temp;
    Circular_kernel::Circular_arc_point_2 p0,p1;
    is >> l_temp >> p0 >> p1;
    //std::cout << "got here l_temp p0 p1 " << l_temp << " " << p0 << " " << p1 << std::endl;
    return Circular_kernel::Line_arc_2(l_temp,p0,p1);
  }
  else if (type == 'w' || type == 'W' || type == 'l' || type == 'L')
  {
    Circular_kernel::Point_2 p0,p1;
    is >> p0 >> p1;
    return Circular_kernel::Line_arc_2(p0,p1);
  }
  else if (type == 'v' || type == 'V')
  {
    Circular_kernel::Segment_2 seg;
    is >> seg;
    return Circular_kernel::Line_arc_2(seg);
  }
  std::cout << "should never happen Line_arc_2 " << type <<std::endl;
  return Circular_kernel::Line_arc_2(); //should never happen
}
#endif

#if TEST_TRAITS == CIRCULAR_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS
template <class stream>
Circular_kernel::Circular_arc_2 read_arc(char type,stream & is)
{
  if (type == 'b' || type == 'B')
  {
    Circular_kernel::Circle_2 c_temp,c_temp1,c_temp2;
    bool b1,b2;
    is >> c_temp >> c_temp1 >> b1 >> c_temp2 >> b2 ;
    return Circular_kernel::Circular_arc_2(c_temp,c_temp1,b1,c_temp2,b2);
  }
  else if (type == 'c' || type == 'C')
  {
    Circular_kernel::Circle_2 c_temp;
    Circular_kernel::Circular_arc_point_2 p0,p1;
    is >> c_temp >> p0 >> p1;
    return Circular_kernel::Circular_arc_2(c_temp,p0,p1);
  }
  else if (type == 'd' || type == 'D')
  {
    Circular_kernel::Circle_2 c_temp;
    Circular_kernel::Line_2 l_temp1,l_temp2;
    bool b1,b2;
    is >> c_temp >> l_temp1 >> b1 >> l_temp2 >> b2;
    return Circular_kernel::Circular_arc_2(c_temp,l_temp1,b1,l_temp2,b2);
  }
  else if (type == 'e' || type == 'E')
  {
    Circular_kernel::Circular_arc_2 a_temp;
    Circular_kernel::Circle_2 c_temp;
    bool b1,b2;
    is >> a_temp >> b1 >> c_temp >> b2;
    return Circular_kernel::Circular_arc_2(a_temp,b1,c_temp,b2);
  }
  std::cout << "should never happen Circular_arc_2" << std::endl;
  return Circular_kernel::Circular_arc_2(); //should never happen
}
#endif

#if TEST_TRAITS == LINE_ARC_TRAITS

/*! Read a line arc point */
template <>
template <class stream>
bool
Traits_base_test<Traits>::
read_point(stream & is, Point_2 & p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone line arc curve */
template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_xcurve(stream & is, X_monotone_curve_2 & xcv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type))
  {
    xcv=read_line(type,is);
    return true;
  }
  return false;

}

/*! Read a general line arc curve */
template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_curve(stream & is, Curve_2 & cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type))
  {
    cv=read_line(type,is);
    return true;
  }
  return false;
}

#elif TEST_TRAITS == CIRCULAR_ARC_TRAITS

/*! Read a circular arc point */
template <>
template <class stream>
bool
Traits_base_test<Traits>::
read_point(stream & is,
           Point_2 & p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone circular arc curve */
template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_xcurve(stream & is,X_monotone_curve_2 & xcv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_2(type))
  {
    xcv=read_arc(type,is);
    return true;
  }
  return false;
}

/*! Read a general circular curve */
template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_curve(stream & is, Curve_2 & cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'a' || type == 'A') 
  {
    Circular_kernel::Circle_2 c_temp;
    is >> c_temp;
    cv=Circular_kernel::Circular_arc_2(c_temp);
    return true;
  }
  else if (is_deg_2(type))
  {
    cv=read_arc(type,is);
    return true;
  }
  return false;
}

#elif TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

/*! Read a circular-line arc point */
template <>
template <class stream>
bool
Traits_base_test<Traits>::
read_point(stream & is, Point_2 & p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone circular-line arc curve */
template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_xcurve(stream & is,X_monotone_curve_2 & xcv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type))
  {
    xcv = read_line(type,is);
    return true;
  }
  else if (is_deg_2(type))
  {
    xcv = X_monotone_curve_2(read_arc(type,is));
    return true;
  }
  return false;
}

/*! Read a general circular-line curve */
template <>
template <class stream>
bool
Traits_base_test<Traits >::
read_curve(stream & is, Curve_2 & cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'a' || type == 'A') 
  {
    Circular_kernel::Circle_2 c_temp;
    is >> c_temp;
    cv=Curve_2(c_temp);
    return true;
  }
  else if (is_deg_1(type))
  {
    cv = Curve_2(read_line(type,is));
    return true;
  }
  else if (is_deg_2(type))
  {
    cv = read_arc(type,is);
    return true;
  }
  return false;

}

#endif

#endif

#endif
