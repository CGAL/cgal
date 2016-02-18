#ifndef CGAL_IO_BASE_TEST_H
#define CGAL_IO_BASE_TEST_H

template <typename GeomTraits_>
class IO_base_test {
public:
  typedef GeomTraits_                                   Geom_traits;
  typedef typename Geom_traits::Point_2                 Point_2;
  typedef typename Geom_traits::X_monotone_curve_2      X_monotone_curve_2;
  typedef typename Geom_traits::Curve_2                 Curve_2;

#if TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS || \
    TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS || \
    TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS || \
    TEST_GEOM_TRAITS == POLYLINE_GEOM_TRAITS ||\
    TEST_GEOM_TRAITS == NON_CACHING_POLYLINE_GEOM_TRAITS
  // Poly curves needs some testing where Segments and X-monotone segments are
  // required instead of polycurves/x-monotone polycurves.
  template <typename InputStream_>
  bool read_segment(InputStream_& is, Subcurve_2& seg);
  template <typename InputStream_>
  bool read_xsegment(InputStream_& is, X_monotone_subcurve_2& xseg);
#endif

  /*! Constructor */
  IO_base_test(const Geom_traits& traits);

  /*! Destructor */
  virtual ~IO_base_test() {}

  template <typename InputStream_>
  bool read_point(InputStream_& is, Point_2&);

  template <typename InputStream_>
  bool read_xcurve(InputStream_& is, X_monotone_curve_2&);

  template <typename InputStream_>
  bool read_curve(InputStream_& is, Curve_2&);

protected:
  /*! An instance of the traits */
  const Geom_traits& m_geom_traits;
};

/*!
 * Constructor.
 * Accepts test data file name.
 */
template <typename GeomTraits_>
IO_base_test<GeomTraits_>::IO_base_test(const GeomTraits_& geom_traits) :
  m_geom_traits(geom_traits) {}

// Generic implementation
template <typename GeomTraits_>
template <typename InputStream_>
bool IO_base_test<GeomTraits_>::
read_point(InputStream_& is, typename GeomTraits_::Point_2& p)
{
  Basic_number_type x, y;
  is >> x >> y;
  p = typename GeomTraits_::Point_2(x, y);
  return true;
}

template <typename GeomTraits_>
template <typename InputStream_>
bool IO_base_test<GeomTraits_>::
read_xcurve(InputStream_& is, typename GeomTraits_::X_monotone_curve_2& xcv)
{
  typedef GeomTraits_   Geom_traits;
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  CGAL_assertion(!is.fail());
  Point_2 p1(x1, y1);
  Point_2 p2(x2, y2);
  CGAL_assertion(p1 != p2);
  xcv = typename Geom_traits::X_monotone_curve_2(p1, p2);
  return true;
}

template <typename GeomTraits_>
template <typename InputStream_>
bool IO_base_test<GeomTraits_>::
read_curve(InputStream_& is, typename GeomTraits_::Curve_2& cv)
{
  typedef GeomTraits_   Geom_traits;
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Point_2 p1(x1, y1);
  Point_2 p2(x2, y2);
  CGAL_assertion(p1 != p2);
  cv = typename Geom_traits::Curve_2(p1, p2);
  return true;
}

// Specialized implementations

// Linear
#if TEST_GEOM_TRAITS == LINEAR_GEOM_TRAITS

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  is >> xcv;
  return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  is >> cv;
  return true;
}

// Polyline
#elif (TEST_GEOM_TRAITS == POLYLINE_GEOM_TRAITS) || \
  (TEST_GEOM_TRAITS == NON_CACHING_POLYLINE_GEOM_TRAITS)

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  unsigned int num_points;
  is >> num_points;
  std::vector<Point_2> points;
  points.clear();
  for (unsigned int j = 0; j < num_points; ++j) {
    Basic_number_type x, y;
    is >> x >> y;
    Point_2 p(x, y);
    points.push_back(p);
  }
  xcv = m_geom_traits.construct_x_monotone_curve_2_object()(points.begin(),
                                                            points.end());
  return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  unsigned int num_points;
  is >> num_points;
  std::vector<Point_2> points;
  points.clear();
  for (unsigned int j = 0; j < num_points; ++j) {
    Basic_number_type x, y;
    is >> x >> y;
    Point_2 p(x, y);
    points.push_back(p);
  }
  cv = m_geom_traits.construct_curve_2_object()(points.begin(), points.end());
  return true;
}
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_segment(InputStream_& is,
                                                  Subcurve_2& seg)
{
  Basic_number_type x, y;
  is >> x >> y;
  Point_2 p_src(x, y);
  is >> x >> y;
  Point_2 p_tgt(x, y);
  seg = Subcurve_2(p_src, p_tgt);
  return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xsegment(InputStream_& is,
                                                   X_monotone_subcurve_2& xseg)
{
  Basic_number_type x, y;
  is >> x >> y;
  Point_2 p_src(x, y);
  is >> x >> y;
  Point_2 p_tgt(x, y);
  xseg = X_monotone_subcurve_2(p_src, p_tgt);
  return true;
}

//polycurve_conic
#elif TEST_GEOM_TRAITS == POLYCURVE_CONIC_GEOM_TRAITS

/*! Read a point */

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_point(InputStream_& is, Point_2& p)
{
  Rational rat_x, rat_y;
  is >> rat_x >> rat_y;
  Basic_number_type x(rat_x), y(rat_y);
  p = Point_2(x, y);
  return true;
}


/*! */
template <typename InputStream_>
bool read_orientation(InputStream_& is, CGAL::Orientation& orient)
{
  int i_orient;
  is >> i_orient;
  orient = (i_orient > 0) ? CGAL::COUNTERCLOCKWISE :
    (i_orient < 0) ? CGAL::CLOCKWISE : CGAL::COLLINEAR;
  return true;
}

/*! */
template <typename InputStream_>
bool read_app_point(InputStream_& is, Point_2& p)
{
  double x, y;
  is >> x >> y;
  p = Point_2(Algebraic(x), Algebraic(y));
  return true;
}

/*! */
template <typename InputStream_>
bool read_orientation_and_end_points(InputStream_& is,
                                     CGAL::Orientation& orient,
                                     Point_2& source, Point_2& target)
{
  // Read the orientation.
  if (!read_orientation(is, orient)) return false;

  // Read the end points of the arc and create it.
  if (!read_app_point(is, source)) return false;
  if (!read_app_point(is, target)) return false;
  return true;
}

/*! */
template <typename InputStream_, typename Curve>
bool read_general_arc(InputStream_& is, Curve& cv)
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
template <typename InputStream_, typename Curve>
bool read_general_conic(InputStream_& is, Curve& cv)
{
  // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
  Rational r, s, t, u, v, w;
  is >> r >> s >> t >> u >> v >> w;
  // Create a full conic (should work only for ellipses).
  cv = Curve(r, s, t, u, v, w);
  return true;
}

/*! */
template <typename InputStream_, typename Curve>
bool read_general_curve(InputStream_& is, Curve& cv)
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

/*! Read an x-monotone conic poly-curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  // since we are dealing with polycurve, we will make more than 1 conic curves
  // (polycurve compatible) and return the x-monotone-constructed polycurve.

  // to store x-monotoneConic curves i.e in Arr_polyline_traits_2 they are
  // called X_monotone_subcurve_2
  std::vector<X_monotone_subcurve_2> conic_x_monotone_segments;

  Subcurve_2 tmp_cv;

  // Get the arc type:
  char type;
  is >> type;

  //get number of x-monotone conic-arcs.
  unsigned int number_of_curves;
  is >> number_of_curves;

  for (unsigned int i=0; i<number_of_curves; ++i) {
    if ((type == 'a') || (type == 'A')) {
      if (!read_general_curve(is, tmp_cv)) return false;
      X_monotone_subcurve_2 tmp_xcv(tmp_cv);
      conic_x_monotone_segments.push_back(tmp_xcv);
    }
    else if ((type == 'c') || (type == 'C')) {
      if (!read_general_conic(is, tmp_cv)) return false;
      X_monotone_subcurve_2 tmp_xcv(tmp_cv);
      conic_x_monotone_segments.push_back(tmp_xcv);
    }
    else if ((type == 'i') || (type == 'I')) {
      if (!read_general_arc(is, tmp_cv)) return false;
      X_monotone_subcurve_2 tmp_xcv(tmp_cv);
      conic_x_monotone_segments.push_back(tmp_xcv);
    }
    else {
      std::cerr << "Illegal conic type specification: " << type << "."
                << std::endl;
      return false;
    }

  } //for loop

  //construct x-monotone polycurve
  xcv = m_geom_traits.construct_x_monotone_curve_2_object()
    (conic_x_monotone_segments.begin(), conic_x_monotone_segments.end());

  return true;
}

/*! Read a conic poly-curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  // since we are dealing with polycurve, we will make more than 1 conic curves
  // (polycurve compatible) and return the constructed polycurve.

  //to store Conic curves i.e in Arr_polyline_traits_2 they are called Subcurve_2
  std::vector<Subcurve_2> conic_segments;

  Subcurve_2 tmp_cv;

  // Get the arc type:
  char type;
  is >> type;

  //get number of xmonotone-conic arcs.
  unsigned int number_of_curves;
  is >> number_of_curves;

  for (unsigned int i = 0; i < number_of_curves; ++i) {
    if ((type == 'a') || (type == 'A')) {
      if (!read_general_curve(is, tmp_cv)) return false;
      conic_segments.push_back(tmp_cv);
    }
    else if ((type == 'c') || (type == 'C')) {
      if (!read_general_conic(is, tmp_cv)) return false;
      conic_segments.push_back(tmp_cv);
    }
    else if ((type == 'i') || (type == 'I')) {
      if (!read_general_arc(is, tmp_cv)) return false;
      conic_segments.push_back(tmp_cv);
    }

    // Enable these later if needed.
    //else if ((type == 'e') || (type == 'E'))
    //  return read_partial_ellipse(is, cv);
    //else if ((type == 'h') || (type == 'H')) return read_hyperbola(is, cv);
    //else if ((type == 'p') || (type == 'P')) return read_parabola(is, cv);
    //else if ((type == 'f') || (type == 'F')) return read_full_ellipse(is, cv);
    //else if ((type == 's') || (type == 'S')) return read_segment(is, cv);

    // If we reached here, we have an unknown conic type:
    else {
      std::cerr << "Illegal conic type specification: " << type << "."
                << std::endl;
      return false;
    }

  } //for loop

  //construct the polycurve
  cv = m_geom_traits.construct_curve_2_object()(conic_segments.begin(),
                                                conic_segments.end());

  return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_segment(InputStream_& is,
                                                  Subcurve_2& seg)
{
  Subcurve_2 tmp_seg;
  char type;
  is >> type;
  if (!read_general_curve(is, tmp_seg)) return false;
  seg = tmp_seg;
   return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xsegment(InputStream_& is,
                                                   X_monotone_subcurve_2& xseg)
{
  char type;
  is >> type;
  Subcurve_2 tmp_seg;
  if (!read_general_curve(is, tmp_seg)) return false;
  xseg = X_monotone_subcurve_2(tmp_seg);
  return true;
}

#elif TEST_GEOM_TRAITS == POLYCURVE_CIRCULAR_ARC_GEOM_TRAITS

/*! */
template <typename InputStream_>
bool read_orientation(InputStream_& is, CGAL::Orientation& orient)
{
  int i_orient;
  is >> i_orient;
  orient = (i_orient > 0) ? CGAL::COUNTERCLOCKWISE :
    (i_orient < 0) ? CGAL::CLOCKWISE : CGAL::COLLINEAR;
  return true;
}

/*! Read an x-monotone circle segment polycurve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  std::vector<X_monotone_subcurve_2> x_segments;
  char type;
  is >> type;
  unsigned int number_of_segments;
  is >> number_of_segments;
  CGAL::Orientation orientation;
  for (unsigned int i = 0; i < number_of_segments; ++i) {
    if ((type == 'x') || (type == 'X')) {
      Rat_point_2 circle_center;
      Point_2 ps, pt;
      Rat_nt circle_radius;

      int point_x, point_y;
      is >> point_x >> point_y;
      circle_center = Rat_point_2(point_x, point_y);

      is >> circle_radius;

      if (!read_orientation(is, orientation)) return false;

      Circle_2 c = Circle_2(circle_center, circle_radius, orientation);

      is >> point_x >> point_y;
      ps = Point_2(Number_type(point_x, 1), Number_type(point_y, 1));

      is >> point_x >> point_y;
      pt = Point_2(Number_type(point_x, 1), Number_type(point_y, 1));

      X_monotone_subcurve_2 x_seg(c, ps, pt, c.orientation());
      x_segments.push_back(x_seg);
    }
    else {
      std::cerr << "Illegal Circle segment type specification: " << type
                << "." << std::endl;
      return false;
    }
  } //for loop

  //construct x-monotone polycurve
  xcv = m_geom_traits.construct_x_monotone_curve_2_object()(x_segments.begin(),
                                                            x_segments.end());

  return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  std::vector<Subcurve_2> segments;
  char type;
  is >> type;
  unsigned int number_of_segments;
  is >> number_of_segments;
  CGAL::Orientation orientation;
  for (unsigned int i = 0; i<  number_of_segments; ++i) {
    if ((type == 'c') || (type == 'C')) {
      Rat_point_2 circle_center;
      Point_2 ps, pt;
      Rat_nt circle_radius;

      int point_x, point_y;
      is >> point_x >> point_y;
      circle_center = Rat_point_2(point_x, point_y);

      is >> circle_radius;

      if (!read_orientation(is, orientation)) return false;

      is >> point_x >> point_y;
      ps = Point_2(Number_type(point_x, 1), Number_type(point_y, 1));

      is >> point_x >> point_y;
      pt = Point_2(Number_type(point_x, 1), Number_type(point_y, 1));

      Subcurve_2 tmp_seg(circle_center, circle_radius, orientation, ps, pt);
      segments.push_back(tmp_seg);
    }
    else {
      std::cerr << "Illegal Circle segment type specification: " << type
                << "." << std::endl;
      return false;
    }
  } //for loop

  //construct polycurve
  cv = m_geom_traits.construct_curve_2_object()(segments.begin(),
                                                segments.end());

  return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_segment(InputStream_& is,
                                                  Subcurve_2& seg)
{
  //we dont need to check this type as it has already been checked in the
  //IO_test.h
  char type;
  is >> type;

  Rat_point_2 circle_center;
  Point_2 ps, pt;
  Rat_nt circle_radius;

  int point_x, point_y;
  is >> point_x >> point_y;
  circle_center = Rat_point_2(point_x, point_y);

  is >> circle_radius;

  CGAL::Orientation orientation;

  if (!read_orientation(is, orientation)) return false;

  is >> point_x >> point_y;
  ps = Point_2(Number_type(point_x, 1), Number_type(point_y, 1));

  is >> point_x >> point_y;
  pt = Point_2(Number_type(point_x, 1), Number_type(point_y, 1));

  Subcurve_2 tmp_seg(circle_center, circle_radius, orientation, ps, pt);

  seg = tmp_seg;

  return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xsegment(InputStream_& is,
                                                   X_monotone_subcurve_2& xseg)
{
  //we dont need to check this type as it has already been checked in the
  //IO_test.h
  char type;
  is >> type;

  CGAL::Orientation orientation;
  Rat_point_2 circle_center;
  Point_2 ps, pt;
  Rat_nt circle_radius;
  int point_x, point_y;
  is >> point_x >> point_y;
  circle_center = Rat_point_2(point_x, point_y);
  is >> circle_radius;
  if (!read_orientation(is, orientation)) return false;

  Circle_2 c = Circle_2(circle_center, circle_radius, orientation);
  is >> point_x >> point_y;
  ps = Point_2(Number_type(point_x, 1), Number_type(point_y, 1));
  is >> point_x >> point_y;
  pt = Point_2(Number_type(point_x, 1), Number_type(point_y, 1));
  X_monotone_subcurve_2 x_seg(c, ps, pt, c.orientation());
  xseg = x_seg;
  return true;
}

#elif TEST_GEOM_TRAITS == POLYCURVE_BEZIER_GEOM_TRAITS

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_point(InputStream_& is, Point_2& p)
{
  char type;
  is >> type;

  //Read bezier segment
  Subcurve_2 seg;
  if (!read_segment(is, seg)) return false;

  Point_2 point;
  if (type == 'r') {
    Rational rat_t;
    is >> rat_t;
    point = Point_2(seg, rat_t);
  }
  else if (type == 'a') {
    Algebraic alg_t;
    is >> alg_t;
    point = Point_2(seg, alg_t);
  }

  p = point;
  return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_segment(InputStream_& is,
                                                  Subcurve_2& seg)
{
  char type;
  is >> type;
  Bezier_tratis bezier_traits;
  std::vector<Control_point_2> point_vector;
  unsigned int num_control_points;
  is >> num_control_points;
  point_vector.clear();
  for (unsigned int j = 0; j < num_control_points; ++j) {
    int point_x, point_y;
    is >> point_x >> point_y;
    point_vector.push_back(Control_point_2(point_x, point_y));
  }
  //get the non x-monotone bezier segment
  seg = Subcurve_2(point_vector.begin(), point_vector.end());
  return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xsegment(InputStream_& is,
                                                   X_monotone_subcurve_2& xseg)
{
  char type;
  is >> type;
  Bezier_tratis bezier_traits;
  std::vector<Control_point_2> point_vector;
  unsigned int num_control_points;
  is >> num_control_points;
  point_vector.clear();
  for (unsigned int j = 0; j < num_control_points; ++j) {
    int point_x, point_y;
    is >> point_x >> point_y;
    point_vector.push_back(Control_point_2(point_x, point_y));
  }
  //get the non x-monotone bezier segment
  Subcurve_2 seg (point_vector.begin(), point_vector.end());

  //convert it into x-monotone bezier segment.
  std::vector<CGAL::Object> obj_vector;
  bezier_traits.make_x_monotone_2_object()(seg,
                                           std::back_inserter(obj_vector));
  X_monotone_subcurve_2 x_segment =
    CGAL::object_cast<X_monotone_subcurve_2>((obj_vector[0]));

  xseg = x_segment;

  return true;
}


template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  std::vector<X_monotone_subcurve_2> x_segments;
  std::vector<Control_point_2> point_vector;

  Bezier_tratis bezier_traits;

  char type;
  is >> type;

  unsigned int number_of_segments;
  is >> number_of_segments;

  if ((type == 'x') || (type == 'X')) {
    for (unsigned int i=0; i<number_of_segments; ++i) {
      unsigned int num_control_points;
      is >> num_control_points;

      point_vector.clear();

      for (unsigned int j=0; j<num_control_points; ++j) {
        int point_x, point_y;
        is >> point_x >> point_y;
        point_vector.push_back(Control_point_2(point_x, point_y));
      }
      //get the non x-monotone bezier segment
      Subcurve_2 seg(point_vector.begin(), point_vector.end());

      //convert it into x-monotone bezier segment.
      std::vector<CGAL::Object> obj_vector;
      bezier_traits.make_x_monotone_2_object()(seg,
                                               std::back_inserter(obj_vector));
      X_monotone_subcurve_2 x_seg =
        CGAL::object_cast<X_monotone_subcurve_2>((obj_vector[0]));

      x_segments.push_back(x_seg);

    } //for loop (number of segments)
  }
  else {
    std::cerr << "Illegal Bezier segment type specification: " << type << "."
              << std::endl;
    return false;
  }

  //construct x-monotone polycurve
  xcv = m_geom_traits.construct_x_monotone_curve_2_object()(x_segments.begin(),
                                                            x_segments.end());
  return true;
}

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  std::vector<Subcurve_2> segments;
  std::vector<Control_point_2> point_vector;

  Bezier_tratis bezier_traits;

  char type;
  is >> type;

  unsigned int number_of_segments;
  is >> number_of_segments;

  if ((type == 'c') || (type == 'C')) {
    for (unsigned int i=0; i<number_of_segments; ++i) {
      unsigned int num_control_points;
      is >> num_control_points;

      point_vector.clear();

      for (unsigned int j=0; j<num_control_points; ++j) {
        int point_x, point_y;
        is >> point_x >> point_y;
        point_vector.push_back(Control_point_2(point_x, point_y));
      }
      //get the non x-monotone bezier segment
      Subcurve_2 seg (point_vector.begin(), point_vector.end());

      segments.push_back(seg);

    } //for loop(number of segments)
  }
  else {
    std::cerr << "Illegal Bezier segment type specification: " << type << "."
              << std::endl;
    return false;
  }

  //construct x-monotone polycurve
  cv = m_geom_traits.construct_curve_2_object()(segments.begin(),
                                                segments.end());
  return true;
}

// Circle segment
#elif TEST_GEOM_TRAITS == CIRCLE_SEGMENT_GEOM_TRAITS

template <typename InputStream_>
bool read_ort_point(InputStream_& is, Point_2& p)
{
  bool is_rat;
  typename Point_2::CoordNT ort_x, ort_y;
  Number_type alpha, beta, gamma;
  is >> is_rat;
  if (is_rat) {
    is >> alpha;
    ort_x = Point_2::CoordNT(alpha);
  }
  else {
    is >> alpha >> beta >> gamma;
    ort_x = Point_2::CoordNT(alpha,beta,gamma);
  }
  is >> is_rat;
  if (is_rat) {
    is >> alpha;
    ort_y=Point_2::CoordNT(alpha);
  }
  else {
    is >> alpha >> beta >> gamma;
    ort_y = Point_2::CoordNT(alpha,beta,gamma);
  }
  p = Point_2(ort_x, ort_y);
  return true;
}

/*! Read an x-monotone circle segment curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  bool ans = true;
  char type;
  is >> type;
  if ((type == 'z') || (type == 'Z')) {
    Line_2 l;
    Point_2 ps, pt;
    is >> l;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    xcv = X_monotone_curve_2(l, ps, pt);
    return ans;
  }
  else if ((type == 'y') || (type == 'Y')) {
    Rat_point_2 ps, pt;
    is >> ps >> pt;
    xcv = X_monotone_curve_2(ps, pt);
    return true;
  }
  else if ((type == 'x') || (type == 'X')) {
    Circle_2 c;
    Point_2 ps,pt;
    is >> c;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    xcv = X_monotone_curve_2(c, ps, pt, c.orientation());
    return ans;
  }
  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal circle segment type specification: " << type
            << std::endl;
  return false;
}

/*! Read a general circle segment curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  bool ans = true;
  char type;
  is >> type;
  if ((type == 'a') || (type == 'A')) {
    Rat_point_2 ps, pt;
    is >> ps >> pt;
    Segment_2 s(ps, pt);
    cv = Curve_2(s);
    return true;
  }
  else if ((type == 'b') || (type == 'B')) {
    Rat_point_2 ps, pt;
    is >> ps >> pt;
    cv = Curve_2(ps, pt);
    return true;
  }
  else if ((type == 'c') || (type == 'C')) {
    Line_2 l;
    Point_2 ps, pt;
    is >> l;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv = Curve_2(l, ps, pt);
    return ans;
  }
  else if ((type == 'd') || (type == 'D')) {
    Circle_2 c;
    is >> c;
    cv = Curve_2(c);
    return true;
  }
  else if ((type == 'e') || (type == 'E')) {
    Rat_point_2 p;
    Rat_nt r;
    int orient;
    is >> p >> r >> orient;
    cv = Curve_2(p, r, static_cast<CGAL::Orientation>(orient));
    return true;
  }
  else if ((type == 'f') || (type == 'F')) {
    Circle_2 c;
    Point_2 ps, pt;
    is >> c;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv = Curve_2(c, ps, pt);
    return ans;
  }
  else if ((type == 'g') || (type == 'G')) {
    Rat_point_2 p;
    Rat_nt r;
    int orient;
    Point_2 ps, pt;
    is >> p >> r >> orient;
    ans &= read_ort_point(is, ps);
    ans &= read_ort_point(is, pt);
    cv = Curve_2(p,r,static_cast<CGAL::Orientation>(orient), ps, pt);
    return ans;
  }
  else if ((type == 'h') || (type == 'H')) {
    Rat_point_2 ps, pm, pt;
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
#elif TEST_GEOM_TRAITS == CORE_CONIC_GEOM_TRAITS

// conic traits and rational traits use same number
// type CORE:Expr so this code can be shared

/*! Read a point */

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_point(InputStream_& is, Point_2& p)
{
  Rational rat_x, rat_y;
  is >> rat_x >> rat_y;
  Basic_number_type x(rat_x), y(rat_y);
  p = Point_2(x, y);
  return true;
}

/*! */
template <typename InputStream_>
bool read_orientation(InputStream_& is, CGAL::Orientation& orient)
{
  int i_orient;
  is >> i_orient;
  orient = (i_orient > 0) ? CGAL::COUNTERCLOCKWISE :
    (i_orient < 0) ? CGAL::CLOCKWISE : CGAL::COLLINEAR;
  return true;
}

/*! */
template <typename InputStream_>
bool read_app_point(InputStream_& is, Point_2& p)
{
  double x, y;
  is >> x >> y;
  p = Point_2(Algebraic(x), Algebraic(y));
  return true;
}

/*! */
template <typename InputStream_>
bool read_orientation_and_end_points(InputStream_& is,
                                     CGAL::Orientation& orient,
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
template <typename InputStream_>
bool read_ellipse(InputStream_& is, bool& is_circle, Rat_circle& circle,
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
template <typename InputStream_, typename Curve>
bool read_partial_ellipse(InputStream_& is, Curve& cv)
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
template <typename InputStream_, typename Curve>
bool read_full_ellipse(InputStream_& is, Curve& cv)
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
template <typename InputStream_, typename Curve>
bool read_hyperbola(InputStream_& is, Curve& cv)
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
template <typename InputStream_, typename Curve>
bool read_parabola(InputStream_& is, Curve& cv)
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
template <typename InputStream_, typename Curve>
bool read_segment(InputStream_& is, Curve& cv)
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
template <typename InputStream_, typename Curve>
bool read_general_arc(InputStream_& is, Curve& cv)
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
template <typename InputStream_, typename Curve>
bool read_general_curve(InputStream_& is, Curve& cv)
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
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  Curve_2 tmp_cv;
  if (!read_curve(is, tmp_cv)) return false;
  xcv = X_monotone_curve_2(tmp_cv);
  return true;
}

/*! Read a general conic curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if ((type == 'f') || (type == 'F')) return read_full_ellipse(is, cv);
  else if ((type == 's') || (type == 'S')) return read_segment(is, cv);
  else if ((type == 'i') || (type == 'I')) return read_general_arc(is, cv);
  else if ((type == 'c') || (type == 'C')) {
    // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
    Rational r, s, t, u, v, w;
    is >> r >> s >> t >> u >> v >> w;
    // Create a full conic (should work only for ellipses).
    cv = Curve_2(r, s, t, u, v, w);
    return true;
  }
  else if ((type == 'e') || (type == 'E')) return read_partial_ellipse(is, cv);
  else if ((type == 'h') || (type == 'H')) return read_hyperbola(is, cv);
  else if ((type == 'p') || (type == 'P')) return read_parabola(is, cv);
  else if ((type == 'a') || (type == 'A')) return read_general_curve(is, cv);

  // If we reached here, we have an unknown conic type:
  std::cerr << "Illegal conic type specification: " << type << "."
	    << std::endl;
  return false;
}

// Rational arc
#elif TEST_GEOM_TRAITS == RATIONAL_ARC_GEOM_TRAITS
/*! Read a point */

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_point(InputStream_& is, Point_2& p)
{
  Base_geom_traits::Construct_point_2 construct_point_2 =
    m_geom_traits.construct_point_2_object();

  Rational x, y;
  is >> x >> y ;
  p = construct_point_2(x, y);
  return true;
}

template <typename InputStream_>
bool read_rational_to_real(InputStream_& is, Algebraic_real_1& r)
{
  static Base_geom_traits::Algebraic_kernel_d_1 algebraic_kernel;
  Rational rat;
  is >> rat;
  r = algebraic_kernel.construct_algebraic_real_1_object()(rat);
  return true;
}

template <typename InputStream_>
bool read_coefficients(InputStream_& is, Rat_vector& coeffs)
{
  unsigned int num_coeffs;
  Rational rat;
  is >> num_coeffs;
  coeffs.clear();
  for (unsigned int j = 0; j < num_coeffs; ++j) {
    is >> rat;
    coeffs.push_back(rat);
  }
  return true;
}

/*! Read a xcurve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  //curve constructor
  const Base_geom_traits::Construct_x_monotone_curve_2
    ctr_x_monotone_curve =
    m_geom_traits.construct_x_monotone_curve_2_object();

  // Get the arc type:
  Rat_vector p_coeffs, q_coeffs;
  Algebraic_real_1 src, trg;
  int dir = 0;
  char type;
  is >> type;
  if ((type == 'a') || (type == 'A')) {
    //Default constructor
    xcv = X_monotone_curve_2();
    return true;
  }
  else if ((type == 'b') || (type == 'B')) {
    //Constructor of a whole polynomial curve
    if (read_coefficients(is,p_coeffs))
      xcv = ctr_x_monotone_curve(p_coeffs.begin(), p_coeffs.end());
    else
      return false;
    return true;
  }
  else if ((type == 'c') || (type == 'C')) {
    //Constructor of a polynomial ray
    if (!read_coefficients(is, p_coeffs)) return false;
    if (!read_rational_to_real(is, src)) return false;
    is >> dir;
    xcv = ctr_x_monotone_curve(p_coeffs.begin(), p_coeffs.end(),
                               src, (dir == 0 ? false : true));
    return true;
  }
  else if ((type == 'd') || (type == 'D')) {
    //Constructor of a polynomial arc
    if (!read_coefficients(is,p_coeffs)) return false;
    if (!read_rational_to_real(is,src)) return false;
    if (!read_rational_to_real(is,trg)) return false;
    xcv = ctr_x_monotone_curve(p_coeffs.begin(), p_coeffs.end(), src, trg);
    return true;
  }
  else if ((type == 'e') || (type == 'E')) {
    //Constructor of a whole rational function
    if (!read_coefficients(is,p_coeffs)) return false;
    if (!read_coefficients(is,q_coeffs)) return false;
    xcv = ctr_x_monotone_curve(p_coeffs.begin(), p_coeffs.end(),
                               q_coeffs.begin(), q_coeffs.end());
    return true;
  }
  else if ((type == 'f') || (type == 'F')) {
    //Constructor of a ray of a rational function
    if (!read_coefficients(is,p_coeffs)) return false;
    if (!read_coefficients(is,q_coeffs)) return false;
    if (!read_rational_to_real(is,src)) return false;
    is >> dir;
    xcv = ctr_x_monotone_curve(p_coeffs.begin(), p_coeffs.end(),
                               q_coeffs.begin(), q_coeffs.end(),
                               src, (dir == 0 ? false : true));
    return true;
  }
  else if ((type == 'g') || (type == 'G')) {
    //Constructor of a bounded rational arc
    if (!read_coefficients(is, p_coeffs)) return false;
    if (!read_coefficients(is, q_coeffs)) return false;
    if (!read_rational_to_real(is,src)) return false;
    if (!read_rational_to_real(is,trg)) return false;

    xcv = ctr_x_monotone_curve( p_coeffs.begin(), p_coeffs.end(),
                                q_coeffs.begin(), q_coeffs.end(),
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
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  //curve constructor
  const Base_geom_traits::Construct_curve_2  construct_curve_2 =
    m_geom_traits.construct_curve_2_object();

  // Get the arc type:
  Rat_vector p_coeffs, q_coeffs;
  Algebraic_real_1 src, trg;
  int dir = 0;
  char type;
  is >> type;
  if ((type == 'a') || (type == 'A')) {
    //Default constructor
    cv = Curve_2();
    return true;
  }
  else if ((type == 'b') || (type == 'B')) {
    //Constructor of a whole polynomial curve
    if (read_coefficients(is,p_coeffs))
      cv = construct_curve_2(p_coeffs.begin(),p_coeffs.end());
    else
      return false;
    return true;
  }
  else if ((type == 'c') || (type == 'C')) {
    //Constructor of a polynomial ray
    if (!read_coefficients(is,p_coeffs)) return false;
    if (!read_rational_to_real(is,src)) return false;
    is >> dir;
    cv = construct_curve_2(p_coeffs.begin(), p_coeffs.end(), src,
                           (dir == 0 ? false : true));
    return true;
  }
  else if ((type == 'd') || (type == 'D')) {
    //Constructor of a polynomial arc
    if (!read_coefficients(is,p_coeffs)) return false;
    if (!read_rational_to_real(is,src)) return false;
    if (!read_rational_to_real(is,trg)) return false;
    cv = construct_curve_2(p_coeffs.begin(), p_coeffs.end(), src, trg);
    return true;
  }
  else if ((type == 'e') || (type == 'E')) {
    //Constructor of a whole rational function
    if (!read_coefficients(is,p_coeffs)) return false;
    if (!read_coefficients(is,q_coeffs)) return false;
    cv = construct_curve_2(p_coeffs.begin(), p_coeffs.end(),
                           q_coeffs.begin(), q_coeffs.end());
    return true;
  }
  else if ((type == 'f') || (type == 'F')) {
    //Constructor of a ray of a rational function
    if (!read_coefficients(is, p_coeffs)) return false;
    if (!read_coefficients(is, q_coeffs)) return false;
    if (!read_rational_to_real(is, src)) return false;
    is >> dir;
    cv =construct_curve_2(p_coeffs.begin(), p_coeffs.end(),
                          q_coeffs.begin(), q_coeffs.end(),
                          src, (dir == 0 ? false : true));
    return true;
  }
  else if ((type == 'g') || (type == 'G')) {
    //Constructor of a bounded rational arc
    if (!read_coefficients(is, p_coeffs)) return false;
    if (!read_coefficients(is, q_coeffs)) return false;
    if (!read_rational_to_real(is,src)) return false;
    if (!read_rational_to_real(is,trg)) return false;
    cv = construct_curve_2(p_coeffs.begin(), p_coeffs.end(),
                           q_coeffs.begin(), q_coeffs.end(),
                           src, trg);
    return true;
  }
  // If we reached here, we have an unknown rational arc type:
  std::cerr << "Illegal rational arc type specification: " << type << "."
            << std::endl;
  return (false);
}

// Bezier
#elif TEST_GEOM_TRAITS == BEZIER_GEOM_TRAITS

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_point(InputStream_& is, Point_2& p)
{
  Rational rat_x, rat_y;
  is >> rat_x >> rat_y;
  p = Point_2(rat_x, rat_y);
  return true;
}

/*! Read an x-monotone bezier curve */

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  std::cout << std::endl;
  std::list<CGAL::Object> x_objs;
  std::list<CGAL::Object>::const_iterator xoit;
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
  Base_geom_traits::Make_x_monotone_2 make_x_monotone =
    this->m_geom_traits.make_x_monotone_2_object();
  make_x_monotone(tmp_cv, std::front_inserter (x_objs));
  xoit = x_objs.begin();
  size_t id(0);
  if (!is.eof()) is >> id;
  std::advance(xoit, id);
  if (CGAL::assign(xcv, *xoit))
    return true;
  return false;
}

/*! Read a general bezier curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  is >> cv;
  return true;
}

// Algebraic
#elif TEST_GEOM_TRAITS == ALGEBRAIC_GEOM_TRAITS

#include <CGAL/IO/io.h>

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_point(InputStream_& is, Point_2& p) {
  Base_geom_traits traits;
  Base_geom_traits::Construct_point_2 construct_point_2 =
    traits.construct_point_2_object();
  char type;
  is >> type;
  switch (type) {
   case 'i': {
     int x = 0, y = 0;
     is >> x >> y;
     p = construct_point_2(x, y);
     break;
    }
   case 'b': {
     Base_geom_traits::Bound x,y;
     is >> x >> y;
     p=construct_point_2(x, y);
     break;
    }
   case 'c': {
     Base_geom_traits::Coefficient x, y;
     is >> x >> y;
     p = construct_point_2(x, y);
     break;
    }
   case 'a': {
     Base_geom_traits::Algebraic_real_1 x, y;
     is >> x >> y;
     p = construct_point_2(x, y);
     break;
    }
   case 's': {
     Base_geom_traits::Algebraic_real_1 x;
     is >> x;
     Base_geom_traits::X_monotone_curve_2 xcv;
     CGAL::swallow(is, '(');
     CGAL_assertion_code(bool check = )
     read_xcurve(is, xcv);
     CGAL_assertion(check);

     CGAL::swallow(is, ')');
     p = construct_point_2(x, xcv);
     break;
    }
   case 'g': {
     Base_geom_traits::Algebraic_real_1 x;
     is >> x;
     Base_geom_traits::Curve_2 c;
     CGAL_assertion_code(bool check = )
     read_curve(is, c);
     CGAL_assertion(check);
     int arcno = 0;
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

/*! Read a xcurve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  Base_geom_traits traits;
  Base_geom_traits::Construct_x_monotone_segment_2 construct_segment_2 =
    traits.construct_x_monotone_segment_2_object();
  char type;
  is >> type;
  switch (type) {
   case '1': {
     Curve_2 cv;
     Point_2 end_left,end_right;
     CGAL_assertion_code(bool check = )
     read_curve(is,cv);
     CGAL_assertion(check);
     CGAL::swallow(is, '(');
     CGAL_assertion_code(check = )
     read_point(is,end_left);
     CGAL_assertion(check);
     CGAL::swallow(is, ')');
     CGAL::swallow(is, '(');
     CGAL_assertion_code(check = )
     read_point(is,end_right);
     CGAL_assertion(check);
     CGAL::swallow(is, ')');
     std::vector<Base_geom_traits::X_monotone_curve_2> xcvs;
     construct_segment_2(cv, end_left, end_right, std::back_inserter(xcvs));
     CGAL_assertion(xcvs.size() == 1);
     xcv = xcvs[0];
     break;
    }
   case '2': {
     Curve_2 cv;
     Point_2 p;
     CGAL_assertion_code(bool check = )
     read_curve(is, cv);
     CGAL_assertion(check);
     CGAL::swallow(is, '(');
     CGAL_assertion_code(check = )
     read_point(is, p);
     CGAL_assertion(check);
     CGAL::swallow(is, ')');
     std::string site_of_p_string;
     Base_geom_traits::Site_of_point site_of_p;
     is >> site_of_p_string;
     if (site_of_p_string == "MIN_ENDPOINT") {
       site_of_p = Base_geom_traits::MIN_ENDPOINT;
     }
     else if (site_of_p_string == "MAX_ENDPOINT") {
       site_of_p = Base_geom_traits::MAX_ENDPOINT;
     }
     else {
       CGAL_assertion(site_of_p_string == "POINT_IN_INTERIOR");
       site_of_p = Base_geom_traits::POINT_IN_INTERIOR;
     }
     std::vector<Base_geom_traits::X_monotone_curve_2> xcvs;
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

/*! Read a curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  Base_geom_traits traits;
  Base_geom_traits::Polynomial_2 p;
  Base_geom_traits::Construct_curve_2 construct_curve_2 =
    traits.construct_curve_2_object();
  is >> p;
  cv = construct_curve_2(p);
  return true;
}

// Spherical arc
#elif TEST_GEOM_TRAITS == GEODESIC_ARC_ON_SPHERE_GEOM_TRAITS

/*! Read a point */

template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_point(InputStream_& is, Point_2& p)
{
  Basic_number_type x, y, z;
  is >> x >> y >> z;
  p = Point_2(x, y, z);
  return true;
}

/*! Read a xcurve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  Point_2 p1, p2;
  read_point(is, p1);
  read_point(is, p2);
  CGAL_assertion(p1 != p2);

  unsigned int flag;
  is >> flag;
  if (flag == 1) {
    X_monotone_curve_2::Direction_3 normal;
    is >> normal;
    xcv = X_monotone_curve_2(p1, p2, normal);
  }
  else
    xcv = X_monotone_curve_2(p1, p2);
  return true;
}

/*! Read a curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  Point_2 p1, p2;
  read_point(is, p1);
  read_point(is, p2);
  CGAL_assertion(p1 != p2);
  unsigned int flag;
  is >> flag;
  if (flag == 1) {
    X_monotone_curve_2::Direction_3 normal;
    is >> normal;
    cv = Curve_2(p1, p2, normal);
  }
  else
    cv = Curve_2(p1, p2);
  return true;
}

// circular line arc
#elif TEST_GEOM_TRAITS == LINE_ARC_GEOM_TRAITS || \
  TEST_GEOM_TRAITS == CIRCULAR_ARC_GEOM_TRAITS || \
  TEST_GEOM_TRAITS == CIRCULAR_LINE_ARC_GEOM_TRAITS

/*! Read an arc point */
template <typename Base_geom_traits_T, typename InputStream_>
bool read_arc_point(InputStream_& is, typename Base_geom_traits_T::Point_2& p)
{
  Basic_number_type x, y;
  is >> x >> y;
  Circular_kernel::Point_2 lp(x, y);
  p = typename Base_geom_traits_T::Point_2(lp);
  return true;
}

bool is_deg_1(char c)
{
  return ((c == 'z') || (c == 'Z') || (c == 'y') || (c == 'Y') ||
          (c == 'x') || (c == 'X') || (c == 'w') || (c == 'W') ||
          (c == 'v') || (c == 'V') || (c == 'l') || (c == 'L'));
}

bool is_deg_2(char c)
{
  return ((c == 'b') || (c == 'B') || (c == 'c') || (c == 'C') ||
          (c == 'd') || (c == 'D') || (c == 'e') || (c == 'E'));
}

#if TEST_GEOM_TRAITS == LINE_ARC_GEOM_TRAITS || \
  TEST_GEOM_TRAITS == CIRCULAR_LINE_ARC_GEOM_TRAITS

template <typename InputStream_>
Circular_kernel::Line_arc_2 read_line(char type, InputStream_& is)
{
  if ((type == 'z') || (type == 'Z')) {
    Circular_kernel::Line_2 l_temp;
    Circular_kernel::Circle_2 c_temp1, c_temp2;
    bool b1, b2;
    is >> l_temp >> c_temp1 >> b1 >> c_temp2 >> b2;
    return Circular_kernel::Line_arc_2(l_temp,c_temp1,b1,c_temp2,b2);
  }
  else if ((type == 'y') || (type == 'Y')) {
    Circular_kernel::Line_2 l_temp, l_temp1, l_temp2;
    is >> l_temp >> l_temp1 >> l_temp2;
    return Circular_kernel::Line_arc_2(l_temp,l_temp1,l_temp2);
  }
  else if ((type == 'x') || (type == 'X')) {
    Circular_kernel::Line_2 l_temp;
    Circular_kernel::Circular_arc_point_2 p0, p1;
    is >> l_temp >> p0 >> p1;
    return Circular_kernel::Line_arc_2(l_temp, p0, p1);
  }
  else if ((type == 'w') || (type == 'W') || (type == 'l') || (type == 'L')) {
    Circular_kernel::Point_2 p0,p1;
    is >> p0 >> p1;
    return Circular_kernel::Line_arc_2(p0, p1);
  }
  else if ((type == 'v') || (type == 'V')) {
    Circular_kernel::Segment_2 seg;
    is >> seg;
    return Circular_kernel::Line_arc_2(seg);
  }
  std::cout << "should never happen Line_arc_2 " << type <<std::endl;
  return Circular_kernel::Line_arc_2(); //should never happen
}

#endif

#if TEST_GEOM_TRAITS == CIRCULAR_ARC_GEOM_TRAITS || \
  TEST_GEOM_TRAITS == CIRCULAR_LINE_ARC_GEOM_TRAITS

template <typename InputStream_>
Circular_kernel::Circular_arc_2 read_arc(char type, InputStream_& is)
{
  if ((type == 'b') || (type == 'B')) {
    Circular_kernel::Circle_2 circle, circle1, circle2;
    bool b1, b2;
    is >> circle >> circle1 >> b1 >> circle2 >> b2;

    return Circular_kernel::Circular_arc_2(circle, circle1, b1, circle2, b2);
  }
  else if ((type == 'c') || (type == 'C')) {
    Circular_kernel::Circle_2 circle;
    Circular_kernel::Circular_arc_point_2 p0, p1;
    is >> circle >> p0 >> p1;
    return Circular_kernel::Circular_arc_2(circle, p0, p1);
  }
  else if ((type == 'd') || (type == 'D')) {
    Circular_kernel::Circle_2 circle;
    Circular_kernel::Line_2 line1, line2;
    bool b1, b2;
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




#if TEST_GEOM_TRAITS == LINE_ARC_GEOM_TRAITS

/*! Read a line arc point */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_point(InputStream_& is, Point_2& p)
{
  typedef InputStream_  Input_stream;
  return read_arc_point<Base_geom_traits, Input_stream>(is, p);
}

/*! Read an x-monotone line arc curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type)) {
    xcv = read_line(type, is);
    return true;
  }
  return false;
}

/*! Read a general line arc curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if (is_deg_1(type)) {
    cv = read_line(type, is);
    return true;
  }
  return false;
}

#endif

#if TEST_GEOM_TRAITS == CIRCULAR_ARC_GEOM_TRAITS

/*! Read a circular arc point */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_point(InputStream_& is, Point_2& p)
{
  typedef InputStream_  Input_stream;
  return read_arc_point<Base_geom_traits, Input_stream>(is, p);
}

/*! Read an x-monotone circular arc curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
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
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if ((type == 'a') || (type == 'A')) {
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

#if TEST_GEOM_TRAITS == CIRCULAR_LINE_ARC_GEOM_TRAITS

/*! Read a circular-line arc point */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_point(InputStream_& is, Point_2& p)
{
  typedef InputStream_  Input_stream;
  return read_arc_point<Base_geom_traits, Input_stream>(is, p);
}

/*! Read an x-monotone circular-line arc curve */
template <>
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_xcurve(InputStream_& is,
                                                 X_monotone_curve_2& xcv)
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
template <typename InputStream_>
bool IO_base_test<Base_geom_traits>::read_curve(InputStream_& is, Curve_2& cv)
{
  // Get the arc type:
  char type;
  is >> type;
  if ((type == 'a') || (type == 'A')) {
    Circular_kernel::Circle_2 circle;
    is >> circle;
    cv = Curve_2(circle);
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
