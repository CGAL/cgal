#ifndef CGAL_IO_TEST_H
#define CGAL_IO_TEST_H

template <typename T_Traits>
class IO_test {
public:
  typedef T_Traits                                      Traits;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Curve_2                      Curve_2;

  template <class stream>
  bool read_point(stream& is, Point_2&);

  template <class stream>
  bool read_xcurve(stream& is, X_monotone_curve_2&);

  template <class stream>
  bool read_curve(stream& is, Curve_2&);

protected:
  /*! An instance of the traits */
  Traits m_traits;
};

// Generic implementation
template <typename T_Traits>
template <typename stream>
bool IO_test<T_Traits>::
read_point(stream& is, typename T_Traits::Point_2& p)
{
  Basic_number_type x, y;
  is >> x >> y;
  p = typename T_Traits::Point_2(x, y);
  return true;
}

template <typename T_Traits>
template <typename stream>
bool IO_test<T_Traits>::
read_xcurve(stream& is, typename T_Traits::X_monotone_curve_2& xcv)
{
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Point_2 p1(x1, y1);
  Point_2 p2(x2, y2);
  CGAL_assertion(p1 != p2);
  xcv = typename T_Traits::X_monotone_curve_2(p1, p2);
  return true;
}

template <typename T_Traits>
template <typename stream>
bool
IO_test<T_Traits>::read_curve(stream& is, typename T_Traits::Curve_2& cv)
{
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Point_2 p1(x1, y1);
  Point_2 p2(x2, y2);
  CGAL_assertion(p1 != p2);
  cv = typename T_Traits::Curve_2(p1, p2);
  return true;
}

// Specialized implementations

// Linear
#if TEST_TRAITS == LINEAR_TRAITS

template <>
template <typename stream>
bool
IO_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  is >> xcv;
  return true;
}

template <>
template <typename stream>
bool
IO_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  is >> cv;
  return true;
}

// Polyline
#elif TEST_TRAITS == POLYLINE_TRAITS || \
  TEST_TRAITS == NON_CACHING_POLYLINE_TRAITS

template <>
template <typename stream>
bool
IO_test<Traits>::
read_xcurve(stream& is, Traits::X_monotone_curve_2& xcv)
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
template <typename stream>
bool
IO_test<Traits>::read_curve(stream& is, Traits::Curve_2& cv)
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
  cv = CGAL::Arr_polyline_traits_2<Segment_traits>::Curve_2(points.begin(),
                                                            points.end());
  return true;
}

// Circle segment
#elif TEST_TRAITS == CIRCLE_SEGMENT_TRAITS

template <typename stream>
bool read_ort_point(stream& is, Point_2& p)
{
  bool is_rat;
  typename Point_2::CoordNT ort_x, ort_y;
  Number_type alpha,beta,gamma;
  is >> is_rat;
  if (is_rat) {
    is >> alpha;
    ort_x=Point_2::CoordNT(alpha);
  }
  else {
    is >> alpha >> beta >> gamma;
    ort_x=Point_2::CoordNT(alpha,beta,gamma);
  }
  is >> is_rat;
  if (is_rat) {
    is >> alpha;
    ort_y=Point_2::CoordNT(alpha);
  }
  else {
    is >> alpha >> beta >> gamma;
    ort_y=Point_2::CoordNT(alpha,beta,gamma);
  }
  p = Point_2(ort_x, ort_y);
  return true;
}

/*! Read an x-monotone circle segment curve */
template <>
template <typename stream>
bool
IO_test<Traits>::read_xcurve(stream& is,X_monotone_curve_2& xcv)
{
  bool ans=true;
  char type;
  is >> type;
  if (type == 'z' || type == 'Z') {
    Line_2 l;
    Point_2 ps,pt;
    is >> l;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    xcv=X_monotone_curve_2(l, ps, pt);
    return ans;
  }
  else if (type == 'y' || type == 'Y') {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    xcv=X_monotone_curve_2(ps,pt);
    return true;
  }
  else if (type == 'x' || type == 'X') {
    Circle_2 c;
    Point_2 ps,pt;
    is >> c;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    xcv=X_monotone_curve_2(c, ps, pt, c.orientation());
    return ans;
  }
  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal circle segment type specification: " << type
            << std::endl;
  return false;
}

/*! Read a general circle segment curve */
template <>
template <typename stream>
bool IO_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  std::cout << "Y1" << std::endl;
  bool ans = true;
  char type;
  is >> type;
  if (type == 'a' || type == 'A') {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    Segment_2 s(ps,pt);
    cv = Curve_2(s);
    return true;
  }
  else if (type == 'b' || type == 'B') {
    Rat_point_2 ps,pt;
    is >> ps >> pt;
    cv = Curve_2(ps,pt);
    return true;
  }
  else if (type == 'c' || type == 'C') {
    Line_2 l;
    Point_2 ps, pt;
    is >> l;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv = Curve_2(l, ps, pt);
    return ans;
  }
  else if (type == 'd' || type == 'D') {
    std::cout << "X1" << std::endl;
    Circle_2 c;
    is >> c;
    std::cout << "X2" << std::endl;
    cv = Curve_2(c);
    std::cout << "X3" << std::endl;
    return true;
  }
  else if (type == 'e' || type == 'E') {
    Rat_point_2 p;
    Rat_nt r;
    int orient;
    is >> p >> r >> orient;
    cv = Curve_2(p,r,static_cast<CGAL::Orientation>(orient));
    return true;
  }
  else if (type == 'f' || type == 'F') {
    Circle_2 c;
    Point_2 ps,pt;
    is >> c;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv = Curve_2(c,ps,pt);
    return ans;
  }
  else if (type == 'g' || type == 'G') {
    Rat_point_2 p;
    Rat_nt r;
    int orient;
    Point_2 ps,pt;
    is >> p >> r >> orient;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv = Curve_2(p,r,static_cast<CGAL::Orientation>(orient), ps, pt);
    return ans;
  }
  else if (type == 'h' || type == 'H') {
    Rat_point_2 ps,pm,pt;
    is >> ps >> pm >> pt;
    cv = Curve_2(ps, pm, pt);
    return true;
  }
  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal circle segment type specification: " << type
            << std::endl;
  return false;
}

// Conic
#elif TEST_TRAITS == CORE_CONIC_TRAITS 

// conic traits and rational traits use same number 
// type CORE:Expr so this code can be shared

/*! Read a point */

template <>
template <typename stream>
bool
IO_test<Traits>::read_point(stream& is, Point_2& p)
{
  Rational rat_x,rat_y;
  is >> rat_x >> rat_y;
  Basic_number_type x(rat_x), y(rat_y);
  p = Point_2(x, y);
  return true;
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
bool read_app_point(stream& is, Point_2& p)
{
  double x, y;
  is >> x >> y;
  p = Point_2(Algebraic(x), Algebraic(y));
  return true;
}

/*! */
template <typename stream>
bool read_orientation_and_end_points(stream& is, CGAL::Orientation& orient,
                                     Point_2& source, Point_2& target)
{
  // Read the orientation.
  if (!read_orientation(is, orient)) return false;

  // Read the end points of the arc and create it.
  if (!read_app_point(is, source)) return false;
  if (!read_app_point(is, target)) return false;
  return true;
}

/*! Read a circle or an ellipse */
template <typename stream>
bool read_ellipse(stream& is, bool& is_circle, Rat_circle& circle,
                  Rational& r, Rational& s,
                  Rational& t, Rational& u, Rational& v, Rational& w)
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
    circle = Rat_circle(Rat_point(x0, y0), a);
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
template <typename stream, typename Curve>
bool read_partial_ellipse(stream& is, Curve& cv)
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
  cv = (is_circle) ? Curve(circle, orient, source, target) :
    Curve(r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! */
template <typename stream, typename Curve>
bool read_full_ellipse(stream& is, Curve& cv)
{
  bool is_circle;               // Is this a circle.
  Rat_circle circle;
  Rational r, s, t, u, v, w;
  if (!read_ellipse(is, is_circle, circle, r, s, t, u, v, w))
    return false;

  // Create a full ellipse (or circle).
  cv = (is_circle) ? Curve(circle) : Curve(r, s, t, u, v, w);
  return true;
}

/*! Read a hyperbola */
template <typename stream, typename Curve>
bool read_hyperbola(stream& is, Curve& cv)
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
  cv = Curve(r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! Read a hyperbola */
template <typename stream, typename Curve>
bool read_parabola(stream& is, Curve& cv)
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
  cv = Curve(r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! */
template <typename stream, typename Curve>
bool read_segment(stream& is, Curve& cv)
{

  // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
  Rational x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  
  Rat_point source(x1, y1);
  Rat_point target(x2, y2);
  Rat_segment segment(source, target);

  // Create the segment.
  cv = Curve(segment);
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
  cv = Curve(r, s, t, u, v, w, orient,
              app_source, r1, s1, t1, u1, v1, w1,
              app_target, r2, s2, t2, u2, v2, w2);
  return true;
}

/*! */
template <typename stream, typename Curve>
bool read_general_curve(stream& is, Curve& cv)
{
  Rational r, s, t, u, v, w;                // The conic coefficients.
  // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
  is >> r >> s >> t >> u >> v >> w;
  CGAL::Orientation orient;
  Point_2 source, target;
  if (!read_orientation_and_end_points(is, orient, source, target))
    return false;

  // Create the conic (or circular) arc.
  cv = Curve(r, s, t, u, v, w, orient, source, target);
  return true;
}

/*! Read an x-monotone conic curve */
template <>
template <typename stream>
bool
IO_test<Traits>::
read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  Curve_2 tmp_cv;
  if (!read_curve(is,tmp_cv))
    return false;
  xcv = X_monotone_curve_2(tmp_cv);
  return true;
}

/*! Read a general conic curve */
template <>
template <typename stream>
bool
IO_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'f' || type == 'F') 
    return read_full_ellipse(is, cv);
  else if (type == 's' || type == 'S') 
    return read_segment(is, cv);
  else if (type == 'i' || type == 'I') 
    return read_general_arc(is, cv);
  else if (type == 'c' || type == 'C') {
    // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
    Rational r, s, t, u, v, w;
    is >> r >> s >> t >> u >> v >> w;
    // Create a full conic (should work only for ellipses).
    cv = Curve_2(r, s, t, u, v, w);
    return true;
  }
  else if (type == 'e' || type == 'E') 
    return read_partial_ellipse(is, cv);
  else if (type == 'h' || type == 'H') 
    return read_hyperbola(is, cv);
  else if (type == 'p' || type == 'P') 
    return read_parabola(is, cv);
  else if (type == 'a' || type == 'A') 
    return read_general_curve(is, cv);

  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal conic type specification: " << type << "."
	    << std::endl;
  return false;
}

// Rational arc
#elif TEST_TRAITS == RATIONAL_ARC_TRAITS
/*! Read a point */

template <>
template <typename stream>
bool
IO_test<Traits>::read_point(stream& is, Point_2& p)
{
  Traits::Construct_point_2 construct_point_2 =
    m_traits.construct_point_2_object();

  Rational x, y;
  is >> x >> y ;
  p = construct_point_2(x, y);
  return true;
}

template <typename stream>
bool read_rational_to_real(stream& is, Algebraic_real_1& r)
{
  static Traits::Algebraic_kernel_d_1 algebraic_kernel;
  Rational rat;
  is >> rat;
  r = algebraic_kernel.construct_algebraic_real_1_object()(rat);
  return true;
}

template <typename stream>
bool read_coefficients(stream& is, Rat_vector& coeffs)
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

/*! Read a xcurve */
template <>
template <typename stream>
bool
IO_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{  
  //curve constructor  
  const Traits::Construct_x_monotone_curve_2
    ctr_x_monotone_curve =
    m_traits.construct_x_monotone_curve_2_object();
  
  // Get the arc type:
  Rat_vector p_coeffs, q_coeffs;
  Algebraic_real_1 src, trg;
  int dir = 0;
  char type;
  is >> type;
  if (type == 'a' || type == 'A') {
    //Default constructor
    xcv = X_monotone_curve_2();
    return true;
  }
  else if (type == 'b' || type == 'B') {
    //Constructor of a whole polynomial curve
    if (read_coefficients(is,p_coeffs))
      xcv = ctr_x_monotone_curve(p_coeffs.begin(),p_coeffs.end());
    else
      return false;
    return true;
  }
  else if (type == 'c' || type == 'C') {
    //Constructor of a polynomial ray
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_rational_to_real(is,src))
      return false;    
    is >> dir;    
    xcv = ctr_x_monotone_curve(p_coeffs.begin(), p_coeffs.end(),
                               src, (dir == 0 ? false : true));
    return true;
  }
  else if (type == 'd' || type == 'D') {
    //Constructor of a polynomial arc
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_rational_to_real(is,src))
      return false;    
    if (!read_rational_to_real(is,trg))
      return false;    
    xcv = ctr_x_monotone_curve(p_coeffs.begin(),p_coeffs.end(), src, trg);
    return true;
  }
  else if (type == 'e' || type == 'E') {
    //Constructor of a whole rational function
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_coefficients(is,q_coeffs))
      return false;
    xcv = ctr_x_monotone_curve(p_coeffs.begin(), p_coeffs.end(),
                               q_coeffs.begin(),q_coeffs.end());
    return true;
  }
  else if (type == 'f' || type == 'F') {
    //Constructor of a ray of a rational function
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_coefficients(is,q_coeffs))
      return false;
    if (!read_rational_to_real(is,src))
      return false;    
    is >> dir;    
    xcv =ctr_x_monotone_curve(p_coeffs.begin(),p_coeffs.end(), 
                              q_coeffs.begin(),q_coeffs.end(), 
                              src, (dir == 0 ? false : true));
    return true;
  }
  else if (type == 'g' || type == 'G') {
    //Constructor of a bounded rational arc
    if (!read_coefficients(is, p_coeffs))
      return false;
    if (!read_coefficients(is, q_coeffs))
      return false;
    if (!read_rational_to_real(is,src))
      return false;    
    if (!read_rational_to_real(is,trg))
      return false;  

    xcv = ctr_x_monotone_curve( p_coeffs.begin(),p_coeffs.end(), 
                                q_coeffs.begin(),q_coeffs.end(), 
                                src, trg);
    return true;
  }
  // If we reached here, we have an unknown rational arc type:
  std::cerr << "Illegal rational arc type specification: " << type << "."
            << std::endl;
  return (false);
}

/*! Read a curve */
template <>
template <typename stream>
bool IO_test<Traits>::read_curve(stream& is, Curve_2& cv)
{  
  //curve constructor
  const Traits::Construct_curve_2  construct_curve_2 =
    m_traits.construct_curve_2_object();

  // Get the arc type:
  Rat_vector p_coeffs, q_coeffs;
  Algebraic_real_1 src, trg;
  int dir = 0;
  char type;
  is >> type;
  if (type == 'a' || type == 'A') {
    //Default constructor
    cv = Curve_2();
    return true;
  }
  else if (type == 'b' || type == 'B') {
    //Constructor of a whole polynomial curve
    if (read_coefficients(is,p_coeffs))
      cv = construct_curve_2(p_coeffs.begin(),p_coeffs.end());
    else
      return false;
    return true;
  }
  else if (type == 'c' || type == 'C') {
    //Constructor of a polynomial ray
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_rational_to_real(is,src))
      return false;    
    is >> dir;
    cv = construct_curve_2(p_coeffs.begin(), p_coeffs.end(), src,
                           (dir == 0 ? false : true));
    return true;
  }
  else if (type == 'd' || type == 'D') {
    //Constructor of a polynomial arc
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_rational_to_real(is,src))
      return false;    
    if (!read_rational_to_real(is,trg))
      return false;    
    cv = construct_curve_2(p_coeffs.begin(),p_coeffs.end(), src, trg);
    return true;
  }
  else if (type == 'e' || type == 'E') {
    //Constructor of a whole rational function
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_coefficients(is,q_coeffs))
      return false;
    cv = construct_curve_2(p_coeffs.begin(), p_coeffs.end(),
                           q_coeffs.begin(), q_coeffs.end());
    return true;
  }
  else if (type == 'f' || type == 'F') {
    //Constructor of a ray of a rational function
    if (!read_coefficients(is,p_coeffs))
      return false;
    if (!read_coefficients(is,q_coeffs))
      return false;
    if (!read_rational_to_real(is,src))
      return false;    
    is >> dir;  
    cv =construct_curve_2(p_coeffs.begin(),p_coeffs.end(), 
                          q_coeffs.begin(),q_coeffs.end(), 
                          src, (dir == 0 ? false : true));
    return true;
  }
  else if (type == 'g' || type == 'G') {
    //Constructor of a bounded rational arc
    if (!read_coefficients(is, p_coeffs))
      return false;
    if (!read_coefficients(is, q_coeffs))
      return false;
    if (!read_rational_to_real(is,src))
      return false;    
    if (!read_rational_to_real(is,trg))
      return false;        
    cv = construct_curve_2(p_coeffs.begin(),p_coeffs.end(), 
                           q_coeffs.begin(),q_coeffs.end(), 
                           src, trg);
    return true;
  }
  // If we reached here, we have an unknown rational arc type:
  std::cerr << "Illegal rational arc type specification: " << type << "."
            << std::endl;
  return (false);
}

// Bezier
#elif TEST_TRAITS == BEZIER_TRAITS

template <>
template <typename stream>
bool
IO_test<Traits>::read_point(stream& is, Point_2& p)
{
  Rational rat_x,rat_y;
  is >> rat_x >> rat_y;
  p = Point_2(rat_x, rat_y);
  return true;
}

/*! Read an x-monotone bezier curve */

template <>
template <typename stream>
bool
IO_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  std::list<CGAL::Object>                  x_objs;
  std::list<CGAL::Object>::const_iterator  xoit;
  Curve_2 tmp_cv;
  is >> tmp_cv;
  Rational B_psx = Rational(tmp_cv.control_point(0).x());
  Rational B_psy = Rational(tmp_cv.control_point(0).y());
  Rational B_ptx =
    Rational(tmp_cv.control_point(tmp_cv.number_of_control_points()-1).x());
  Rational B_pty =
    Rational(tmp_cv.control_point(tmp_cv.number_of_control_points()-1).y());
  Point_2 B_ps(B_psx, B_psy);
  Point_2 B_pt(B_ptx, B_pty);
  Traits::Make_x_monotone_2 make_x_monotone =
    this->m_traits.make_x_monotone_2_object();
  make_x_monotone (tmp_cv, std::front_inserter (x_objs));
  xoit = x_objs.begin();
  if (CGAL::assign(xcv, *xoit))
    return true;
  return false;
}

/*! Read a general bezier curve */
template <>
template <typename stream>
bool
IO_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  is >> cv;
  return true;
}

// Algebraic
#elif TEST_TRAITS == ALGEBRAIC_TRAITS

#include <CGAL/IO/io.h>

template <>
template <typename stream>
bool
IO_test<Traits>::read_point(stream& is, Point_2& p) {
  Traits traits;
  Traits::Construct_point_2 construct_point_2 =
    traits.construct_point_2_object();
  char type;
  is >> type;
  switch(type) {
   case 'i': {
     int x=0,y=0;
     is >> x >> y;
     p = construct_point_2(x, y);
     break;
    }
   case 'b': {
     Traits::Bound x,y;
     is >> x >> y;
     p=construct_point_2(x, y);
     break;
    }
   case 'c': {
     Traits::Coefficient x, y;
     is >> x >> y;
     p = construct_point_2(x, y);
     break;
    }
   case 'a': {
     Traits::Algebraic_real_1 x, y;
     is >> x >> y;
     p = construct_point_2(x, y);
     break;
    }
   case 's': {
     Traits::Algebraic_real_1 x;
     is >> x;
     Traits::X_monotone_curve_2 xcv;
     CGAL::swallow(is,'(');
     CGAL_assertion_code(bool check=)
       read_xcurve(is, xcv);
     CGAL_assertion(check);
    
     CGAL::swallow(is,')');
     p = construct_point_2(x, xcv);
     break;
    }
   case 'g': {
     Traits::Algebraic_real_1 x;
     is >> x;
     Traits::Curve_2 c;
     CGAL_assertion_code(bool check = )
       read_curve(is,c);
     CGAL_assertion(check);
     int arcno=0;
     is >> arcno;
     p = construct_point_2(x, c, arcno);
     break;
    }
   default: {
     std::cout << "Expected i, b, c, a, s, or g, but got \"" << type << "\""
               << std::endl;
     return false;
   }
  }
  return true;
}

template <>
template <typename stream>
bool IO_test<Traits>::read_xcurve(stream& is, 
                                           Traits::X_monotone_curve_2& xcv)
{
  Traits traits;
  Traits::Construct_x_monotone_segment_2 construct_segment_2 =
    traits.construct_x_monotone_segment_2_object();
  char type;
  is >> type;
  switch(type) {
   case '1': {
     Curve_2 cv;
     Point_2 end_left,end_right;
     CGAL_assertion_code(bool check=)
     read_curve(is,cv);
     CGAL_assertion(check);
     CGAL::swallow(is,'(');
     CGAL_assertion_code(check=)
     read_point(is,end_left);
     CGAL_assertion(check);
     CGAL::swallow(is,')');
     CGAL::swallow(is,'(');
     CGAL_assertion_code(check=)
     read_point(is,end_right);
     CGAL_assertion(check);
     CGAL::swallow(is,')');
     std::vector<Traits::X_monotone_curve_2> xcvs;
     construct_segment_2(cv, end_left, end_right, std::back_inserter(xcvs));
     CGAL_assertion(xcvs.size() == 1);
     xcv = xcvs[0];
     break;
    }
   case '2': {
     Curve_2 cv;
     Point_2 p;
     CGAL_assertion_code(bool check=)
     read_curve(is,cv);
     CGAL_assertion(check);
     CGAL::swallow(is,'(');
     CGAL_assertion_code(check=)
     read_point(is,p);
     CGAL_assertion(check);
     CGAL::swallow(is,')');
     std::string site_of_p_string;
     Traits::Site_of_point site_of_p;
     is >> site_of_p_string;
     if (site_of_p_string=="MIN_ENDPOINT") {
       site_of_p=Traits::MIN_ENDPOINT;
     } else if (site_of_p_string=="MAX_ENDPOINT") {
       site_of_p=Traits::MAX_ENDPOINT;
     } else {
       CGAL_assertion(site_of_p_string=="POINT_IN_INTERIOR");
       site_of_p=Traits::POINT_IN_INTERIOR;
     }
     std::vector<Traits::X_monotone_curve_2> xcvs;
     construct_segment_2(cv, p, site_of_p, std::back_inserter(xcvs));
     CGAL_assertion(xcvs.size() == 1);
     xcv = xcvs[0];
     break;
    }
   default: {
     std::cout << "Expected 1 or 2, but got \"" << type << "\"" << std::endl;
     return false;
   }
  }
  return true;
}

template <>
template <typename stream>
bool IO_test<Traits>::read_curve(stream& is, Curve_2& cv) {
  Traits traits;
  Traits::Polynomial_2 p;
  Traits::Construct_curve_2 construct_curve_2 =
    traits.construct_curve_2_object();
  is >> p;
  cv = construct_curve_2(p);
  return true;
}

// Spherical arc
#elif TEST_TRAITS == SPHERICAL_ARC_TRAITS

/*! Read a point */

template <>
template <typename stream>
bool
IO_test<Traits>::read_point(stream& is, Point_2& p)
{
  Basic_number_type x, y, z;
  is >> x >> y >> z;
  p = Point_2(x, y, z);
  return true;
}

/*! Read a xcurve */
template <>
template <typename stream>
bool
IO_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
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
template <typename stream>
bool
IO_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  Point_2 p1, p2;
  read_point(is, p1);
  read_point(is, p2);
  CGAL_assertion(p1 != p2);
  cv = Curve_2(p1, p2);
  return true;
}

// circular line arc
#elif TEST_TRAITS == LINE_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

/*! Read an arc point */
template <typename T_Traits, typename stream>
bool read_arc_point(stream& is, typename T_Traits::Point_2& p)
{
  Basic_number_type x, y;
  is >> x >> y;
  Circular_kernel::Point_2 lp(x, y);
  p = typename T_Traits::Point_2(lp);
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

template <typename stream>
Circular_kernel::Line_arc_2 read_line(char type, stream& is)
{
  if (type == 'z' || type == 'Z') {
    Circular_kernel::Line_2 l_temp;
    Circular_kernel::Circle_2 c_temp1,c_temp2;
    bool b1,b2;
    is >> l_temp >> c_temp1 >> b1 >> c_temp2 >> b2;
    return Circular_kernel::Line_arc_2(l_temp,c_temp1,b1,c_temp2,b2);
  }
  else if (type == 'y' || type == 'Y') {
    Circular_kernel::Line_2 l_temp,l_temp1,l_temp2;
    is >> l_temp >> l_temp1 >> l_temp2;
    return Circular_kernel::Line_arc_2(l_temp,l_temp1,l_temp2);
  }
  else if (type == 'x' || type == 'X') {
    Circular_kernel::Line_2 l_temp;
    Circular_kernel::Circular_arc_point_2 p0,p1;
    is >> l_temp >> p0 >> p1;
    //std::cout << "got here l_temp p0 p1 " << l_temp << " " << p0 << " " << p1 << std::endl;
    return Circular_kernel::Line_arc_2(l_temp, p0, p1);
  }
  else if (type == 'w' || type == 'W' || type == 'l' || type == 'L') {
    Circular_kernel::Point_2 p0,p1;
    is >> p0 >> p1;
    return Circular_kernel::Line_arc_2(p0,p1);
  }
  else if (type == 'v' || type == 'V') {
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
template <typename stream>
Circular_kernel::Circular_arc_2 read_arc(char type,stream& is)
{
  if (type == 'b' || type == 'B') {
    Circular_kernel::Circle_2 circle, circle1, circle2;
    bool b1, b2;
    is >> circle >> circle1 >> b1 >> circle2 >> b2;

    return Circular_kernel::Circular_arc_2(circle, circle1, b1, circle2, b2);
  }
  else if (type == 'c' || type == 'C') {
    Circular_kernel::Circle_2 circle;
    Circular_kernel::Circular_arc_point_2 p0, p1;
    is >> circle >> p0 >> p1;
    return Circular_kernel::Circular_arc_2(circle, p0, p1);
  }
  else if (type == 'd' || type == 'D') {
    Circular_kernel::Circle_2 circle;
    Circular_kernel::Line_2 line1, line2;
    bool b1,b2;
    is >> circle >> line1 >> b1 >> line2 >> b2;
    return Circular_kernel::Circular_arc_2(circle, line1, b1, line2, b2);
  }
  else
    CGAL_error_msg("Unrecognized constructor. Should never happen"      \
                   "Circular_arc_2");
  //  else if (type == 'e' || type == 'E')
  //  {
  //    Circular_kernel::Circular_arc_2 arc;
  //    Circular_kernel::Circle_2 circle;
  //    bool b1, b2;
  //    is >> arc >> b1 >> circle >> b2;
  //    return Circular_kernel::Circular_arc_2(arc, b1, circle, b2);
  //  }
  return Circular_kernel::Circular_arc_2(); //should never happen
}
#endif

#if TEST_TRAITS == LINE_ARC_TRAITS

/*! Read a line arc point */
template <>
template <typename stream>
bool
IO_test<Traits>::read_point(stream& is, Point_2& p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone line arc curve */
template <>
template <typename stream>
bool
IO_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type)) {
    xcv = read_line(type,is);
    return true;
  }
  return false;
}

/*! Read a general line arc curve */
template <>
template <typename stream>
bool
IO_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type)) {
    cv = read_line(type,is);
    return true;
  }
  return false;
}

#endif

#if TEST_TRAITS == CIRCULAR_ARC_TRAITS

/*! Read a circular arc point */
template <>
template <typename stream>
bool
IO_test<Traits>::read_point(stream& is, Point_2& p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone circular arc curve */
template <>
template <typename stream>
bool
IO_test<Traits>::read_xcurve(stream& is,X_monotone_curve_2& xcv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_2(type)) {
    xcv = read_arc(type, is);
    return true;
  }
  return false;
}

/*! Read a general circular curve */
template <>
template <typename stream>
bool
IO_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'a' || type == 'A') {
    Circular_kernel::Circle_2 circle;
    is >> circle;
    cv = Circular_kernel::Circular_arc_2(circle);
    return true;
  }
  else if (is_deg_2(type)) {
    cv = read_arc(type, is);
    return true;
  }
  return false;
}

#endif

#if TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

/*! Read a circular-line arc point */
template <>
template <typename stream>
bool
IO_test<Traits>::read_point(stream& is, Point_2& p)
{
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone circular-line arc curve */
template <>
template <typename stream>
bool
IO_test<Traits>::read_xcurve(stream& is, X_monotone_curve_2& xcv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type)) {
    xcv = read_line(type,is);
    return true;
  }
  else if (is_deg_2(type)) {
    xcv = X_monotone_curve_2(read_arc(type, is));
    return true;
  }
  return false;
}

/*! Read a general circular-line curve */
template <>
template <typename stream>
bool IO_test<Traits>::read_curve(stream& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (type == 'a' || type == 'A') {
    Circular_kernel::Circle_2 circle;
    is >> circle;
    cv=Curve_2(circle);
    return true;
  }
  else if (is_deg_1(type)) {
    cv = Curve_2(read_line(type, is));
    return true;
  }
  else if (is_deg_2(type)) {
    cv = read_arc(type, is);
    return true;
  }
  return false;

}

#endif

#endif

#endif
