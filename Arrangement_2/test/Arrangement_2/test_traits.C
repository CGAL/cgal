#include <CGAL/basic.h>
#include <CGAL/Arrangement_2.h>

#include <vector>

#include "test_configuration.h"
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
  Traits_test<Traits> test(argc, argv);
  return (test.start()) ? 0 /* SUCCESS */ : -1; /* FAILURE */
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
Traits_test<CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits> >::
read_point(stream & is,
           CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>::
           Point_2 & p)
{
  Basic_number_type x, y;
  is >> x >> y;
  p = CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>::Point_2(x, y);
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
  //
  //     x - x0   2      y - y0   2
  //  ( -------- )  + ( -------- )  = 1
  //       a               b
  //
  Rational a, b, x0, y0;
  is >> a >> b >> x0 >> y0;
    
  Rational a_sq = a*a;
  Rational b_sq = b*b;

  if (a == b) {
    is_circle = true;
    circle = Rat_circle (Traits::Rat_point_2 (x0, y0), a*b);
  }
  else {
    r = b_sq;
    s = a_sq;
    t = 0;
    u = -2*x0*b_sq;
    v = -2*y0*a_sq;
    w = x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq;
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
  if (read_orientation_and_end_points(is, orient, source, target))
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
  //
  //     x - x0   2      y - y0   2
  //  ( -------- )  - ( -------- )  = 1
  //       a               b
  //
  Rational a, b, x0, y0;
  is >> a >> b >> x0 >> y0;

  Rational a_sq = a * a;
  Rational b_sq = b * b;

  
  Rational r = b_sq;
  Rational s= -a_sq;
  Rational t = 0;
  Rational u = -2*x0*b_sq;
  Rational v = 2*y0*a_sq;
  Rational w = x0*x0*b_sq - y0*y0*a_sq - a_sq*b_sq;  

  CGAL::Orientation orient;
  Point_2 source, target;
  if (read_orientation_and_end_points(is, orient, source, target))
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
  //             2
  //  4c*(y - y0) = (x - x0)
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
  if (read_orientation_and_end_points(is, orient, source, target))
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
  if (read_orientation_and_end_points(is, orient, source, target))
    return false;

  // Create the conic (or circular) arc.
  cv = Curve (r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! Read an x-monotone conic curve */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits> >::
read_xcurve(stream & is,
            CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>::
            X_monotone_curve_2 & cv)
{
  // Get the arc type:
  char type;
  is >> type;

  if (0) ;
#if 0
  else if (type == 's' || type == 'S') {
    return read_segment(is, cv);
  }
#endif
#if 0
  else if (type == 'i' || type == 'I') {
    return read_general_arc(is, cv);
  }
#endif
#if 0
  else if (type == 'c' || type == 'C') {
    // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
    Rational r, s, t, u, v, w;
    is >> r >> s >> t >> u >> v >> w;
    // Create a full conic (should work only for ellipses).
    cv = X_monotone_curve_2 (r, s, t, u, v, w);
    return true;
  }
#endif
#if 0
  else if (type == 'e' || type == 'E') {
    return read_partial_ellipse(is, cv);
  }
#endif
#if 0
  else if (type == 'h' || type == 'H') {
    return read_hyperbola(is, cv);
  }
#endif
#if 0
  else if (type == 'p' || type == 'P') {
    return read_parabola(is, cv);
  }
#endif
#if 0
  else if (type == 'a' || type == 'A') {
    return read_general_curve(is, cv);
  }
#endif
  else {
    std::cerr << "Illegal conic type specification: " << type << "."
              << std::endl;
    return false;
  }

  return true;
}

/*! Read a general conic curve */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits> >::
read_curve(stream & is,
           CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>::
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

#endif
