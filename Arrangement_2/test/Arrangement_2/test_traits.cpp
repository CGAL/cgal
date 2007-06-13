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
#include "test_arc_traits.h"

// Arrangement types:
typedef CGAL::Arr_default_dcel<Traits>                  Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel>               Arr;

// Traits types:
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef Traits::Curve_2                                 Curve_2;

int main (int argc, char * argv[])
{
  CGAL::set_error_behaviour(CGAL::CONTINUE);
  CGAL::set_warning_behaviour(CGAL::CONTINUE);
  prev_error_handler = CGAL::set_error_handler(failure_handler);
  prev_warning_handler = CGAL::set_warning_handler(failure_handler);
  Traits_test<Traits> test(argc, argv);
  bool rc;
  rc = test.start();
  CGAL::set_error_handler(prev_error_handler);
  CGAL::set_warning_handler(prev_warning_handler);
  return (rc) ? 0 : -1;
}

#if TEST_TRAITS == POLYLINE_TRAITS || TEST_TRAITS == NON_CACHING_POLYLINE_TRAITS

template <>
template <class stream>
bool
Traits_test<CGAL::Arr_polyline_traits_2<Segment_traits> >::
read_xcurve(stream & is,
            CGAL::
            Arr_polyline_traits_2<Segment_traits>::X_monotone_curve_2 & cv)
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
    Arr_polyline_traits_2<Segment_traits>::X_monotone_curve_2(points.begin(),
                                                              points.end());
  return true;
}

template <>
template <class stream>
bool
Traits_test<CGAL::Arr_polyline_traits_2<Segment_traits> >::
read_curve(stream & is,
           CGAL::Arr_polyline_traits_2<Segment_traits>::Curve_2 & cv)
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

#elif TEST_TRAITS == CORE_CONIC_TRAITS

template <>
template <class stream>
bool
Traits_test<Traits >::
read_point(stream & is,
           Traits::
           Point_2 & p)
{
  Basic_number_type x, y;
  is >> x >> y;
  //std::cout << "exec x = " << x << std::endl << "y = " << y << std::endl;
  p = Point_2(x, y);
  return true;
}

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
  //std::cout << "app x = " << x << " y = " << y << std::endl; 
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

  //Rational a_sq = a*a; ???!!!
  //Rational b_sq = b*b; ???!!!

  if (a == b) {
    is_circle = true;
    circle = Rat_circle (Traits::Rat_point_2 (x0, y0), a);
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
  //std::cout << "r = " << r  << " s = " << s << " t = " << t << " u = " << u << " v = " << v << " w = " << w <<std::endl;
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

  //Rational a_sq = a * a;
  //Rational b_sq = b * b;

  
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
  typedef Traits::Rat_point_2        Rat_point;
  typedef Traits::Rat_segment_2      Rat_segment;

  // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
  Rational x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  
  //std::cout << "x1 " << x1 << " y1 " << y1 << " x2 " << x2 << " y2 " << y2 << std::endl;
  
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
  //std::cout << r <<" "<< s <<" "<< t <<" "<< u <<" "<< v <<" "<< w <<std::endl;
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
Traits_test<Traits>::
read_xcurve(stream & is,
            Traits::X_monotone_curve_2 & xcv)
{
  Traits::Curve_2 tmp_cv;
  if (!read_curve(is,tmp_cv))
    return false;
  xcv=Traits::X_monotone_curve_2(tmp_cv);

//std::cout <<"debug " <<xcv<<std::endl;

  return true;
}

/*! Read a general conic curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_curve(stream & is,
           Traits::Curve_2 & cv)
{
  // Get the arc type:
  char type;
  is >> type;

  //std::cout << "type = " << type << std::endl;

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

#elif TEST_TRAITS == CIRCLE_TRAITS

/* is it in use ??? to ask efi about Has_boundary_category */

/*! Read an x-monotone conic curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_xcurve(stream & is,Traits::X_monotone_curve_2 & xcv)
{
  is >> xcv;
  return true;
}

/*! Read a general conic curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_curve(stream & is,Traits::Curve_2 & cv)
{
  is >> cv;
  return true;
}

#elif TEST_TRAITS == CIRCLE_SEGMENT_TRAITS

/*! Read an x-monotone conic curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_xcurve(stream & is,Traits::X_monotone_curve_2 & xcv)
{
  char type;
  is >> type;
  if (type == 'z' || type == 'Z')
  {
    Line_2 l;
    Point_2 ps,pt;
    is >> l >> ps >> pt;
    xcv=Traits::X_monotone_curve_2(l,ps,pt);
    return true;
  }
  else if (type == 'y' || type == 'Y')
  {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    xcv=Traits::X_monotone_curve_2(ps,pt);
    return true;
  }
  else if (type == 'x' || type == 'X')
  {
    Circle_2 c;
    Point_2 ps,pt;
    is >> c >> ps >> pt;
/*
    std::cout << "c " << c << std::endl;
    std::cout << "ps " << ps << std::endl;
    std::cout << "pt " << pt << std::endl;
    std::cout << "orient " << orient << std::endl;
*/
    xcv=Traits::X_monotone_curve_2(c,ps,pt,c.orientation());
    return true;
  }
  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal circle segment type specification: " << type << std::endl;
  return false;
}

/*! Read a general conic curve */
template <>
template <class stream>
bool
Traits_test<Traits >::
read_curve(stream & is,Traits::Curve_2 & cv)
{
  char type;
  is >> type;
  if (type == 'a' || type == 'A')
  {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    Segment_2 s(ps,pt);
    cv=Traits::Curve_2(s);
    return true;
  }
  else if (type == 'b' || type == 'B')
  {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    cv=Traits::Curve_2(ps,pt);
    return true;
  }
  else if (type == 'c' || type == 'C')
  {
    Line_2 l;
    Point_2 ps,pt;
    is >> l >> ps >> pt;
    cv=Traits::Curve_2(l,ps,pt);
    return true;
  }
  else if (type == 'd' || type == 'D')
  {
    Circle_2 c;
    is >> c;
    cv=Traits::Curve_2(c);
    return true;
  }
  else if (type == 'e' || type == 'E')
  {
    Rat_point_2 p;
    Rat_nt r;
    int orient;
    is >> p >> r >> orient;
    cv=Traits::Curve_2(p,r,static_cast<CGAL::Orientation>(orient));
    return true;
  }
  else if (type == 'f' || type == 'F')
  {
    Circle_2 c;
    Point_2 ps,pt;
    is >> c >> ps >> pt;
    cv=Traits::Curve_2(c,ps,pt);
    return true;
  }
  else if (type == 'g' || type == 'G')
  {
    Rat_point_2 p;
    Rat_nt r;
    int orient;
    Point_2 ps,pt;
    is >> p >> r >> orient >> ps >> pt;
    cv=Traits::Curve_2(p,r,static_cast<CGAL::Orientation>(orient),ps,pt);
    return true;
  }
  else if (type == 'h' || type == 'H')
  {
    Rat_point_2 ps,pm,pt;
    is >> ps >> pm >> pt;
    cv=Traits::Curve_2(ps,pm,pt);
    return true;
  }
  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal circle segment type specification: " << type << std::endl;
  return false;
}

#elif TEST_TRAITS == BEZIER_TRAITS
/*! Read an x-monotone bezier curve */

Traits traits;

Traits::Make_x_monotone_2  make_x_monotone = traits.make_x_monotone_2_object();

std::list<CGAL::Object>                  x_objs;

std::list<CGAL::Object>::const_iterator  xoit;

template <>
template <class stream>
bool
Traits_test<Traits>::read_xcurve(stream & is,Traits::X_monotone_curve_2 & xcv)
{
  Traits::Curve_2 tmp_cv;

  is >> tmp_cv;

  Algebraic B_psx = Algebraic(tmp_cv.control_point(0).x());
  Algebraic B_psy = Algebraic(tmp_cv.control_point(0).y());
  Algebraic B_ptx = Algebraic(tmp_cv.control_point(tmp_cv.number_of_control_points()-1).x());
  Algebraic B_pty = Algebraic(tmp_cv.control_point(tmp_cv.number_of_control_points()-1).y());

  //std::cout << "B_psx " << B_psx << " B_psy " << B_psy << " B_ptx " << B_ptx << " B_pty " << B_pty << std::endl;

  Traits::Point_2 B_ps(B_psx, B_psy);
  Traits::Point_2 B_pt(B_ptx, B_pty);

  //std::cout << "B_ps " << B_ps << " B_pt " << B_pt << std::endl;

  make_x_monotone (tmp_cv, std::back_inserter (x_objs));
  int idx=0;
  for (xoit = x_objs.begin(); xoit != x_objs.end() && idx==0; ++xoit) {
    // this loop should make only one iteration
    idx++;
    if (CGAL::assign (xcv, *xoit))
      return true;
    else
      return false;
  }
  //xcv = Traits::X_monotone_curve_2(tmp_cv, B_ps, B_pt, Traits::Bezier_cache());
  return true;
}

/*! Read a general bezier curve */
template <>
template <class stream>
bool
Traits_test<Traits >::read_curve(stream & is,Traits::Curve_2 & cv)
{
  is >> cv;
  return true;
}

#endif

#endif
