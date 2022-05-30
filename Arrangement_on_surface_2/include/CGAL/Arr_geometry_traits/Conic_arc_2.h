// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_CONIC_ARC_2_H
#define CGAL_CONIC_ARC_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Header file for the _Conic_arc_2<Int_kernel, Alg_kernel, Nt_traits> class.
 */

#include <ostream>

#include <CGAL/Arr_geometry_traits/Conic_point_2.h>
#include <CGAL/Arr_geometry_traits/Conic_intersections_2.h>
#include <CGAL/Bbox_2.h>

namespace CGAL {

/*! Representation of a conic arc -- a bounded segment that lies on a conic
 * curve, the loci of all points satisfying the equation:
 *   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0
 *
 * The class is templated with three parameters:
 * Rat_kernel A kernel that provides the input objects or coefficients.
 *            Rat_kernel::FT should be an integral or a rational type.
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices, which are algebraic
 *            numbers of degree up to 4 (preferably it is CORE::Expr).
 * Nt_traits A traits class for performing various operations on the integer,
 *           rational and algebraic types.
 */

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class _Conic_arc_2 {
public:
  typedef RatKernel                                        Rat_kernel;
  typedef AlgKernel                                        Alg_kernel;
  typedef NtTraits                                         Nt_traits;

  typedef _Conic_arc_2<Rat_kernel, Alg_kernel, Nt_traits>  Self;

  typedef typename Rat_kernel::FT                          Rational;
  typedef typename Rat_kernel::Point_2                     Rat_point_2;
  typedef typename Rat_kernel::Segment_2                   Rat_segment_2;
  typedef typename Rat_kernel::Circle_2                    Rat_circle_2;

  typedef typename Nt_traits::Integer                      Integer;

  typedef typename Alg_kernel::FT                          Algebraic;
  typedef typename Alg_kernel::Point_2                     Point_2;
  typedef _Conic_point_2<Alg_kernel>                       Conic_point_2;

  /*! \struct
   * For arcs whose base is a hyperbola we store the axis (a*x + b*y + c = 0)
   * which separates the two bracnes of the hyperbola. We also store the side
   * (NEGATIVE or POSITIVE) that the arc occupies.
   * In case of line segments connecting two algebraic endpoints, we use this
   * structure two store the coefficients of the line supporting this segment.
   * In this case we set the side field to be ZERO.
   */
  struct Extra_data {
    Algebraic a;
    Algebraic b;
    Algebraic c;
    Sign side;
  };

  // Bit masks for the m_info field.
  enum {
    IS_VALID = 0,
    IS_FULL_CONIC,
    LAST_INFO,
  };

protected:
  Integer m_r;          //
  Integer m_s;          // The coefficients of the supporting conic curve:
  Integer m_t;          //
  Integer m_u;          //
  Integer m_v;          //   r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
  Integer m_w;          //

  Orientation m_orient;         // The orientation of the conic.
  int m_info;                   // does the arc represent a full conic curve.
  Conic_point_2 m_source;       // the source of the arc (if not a full curve).
  Conic_point_2 m_target;       // the target of the arc (if not a full curve).
  Extra_data* m_extra_data;     // The extra data stored with the arc
                                // (may be nullptr).

public:
  /// \name Construction and destruction functions.
  //@{

  /*! Default constructor.
   */
  _Conic_arc_2() :
    m_r(0), m_s(0), m_t(0), m_u(0), m_v(0), m_w(0),
    m_orient(COLLINEAR),
    m_info(0),
    m_extra_data(nullptr)
  {}

  /*! Copy constructor.
   * \param arc The copied arc.
   */
  _Conic_arc_2(const Self& arc) :
    m_r(arc.m_r), m_s(arc.m_s), m_t(arc.m_t),
    m_u(arc.m_u), m_v(arc.m_v), m_w(arc.m_w),
    m_orient(arc.m_orient),
    m_info(arc.m_info),
    m_source(arc.m_source),
    m_target(arc.m_target)
  {
    m_extra_data = (arc.m_extra_data != nullptr) ?
      new Extra_data(*(arc.m_extra_data)) : nullptr;
  }

  /*! Destructor.
   */
  virtual ~_Conic_arc_2() { if (m_extra_data != nullptr) delete m_extra_data; }

  /*! Assignment operator.
   * \param arc The copied arc.
   */
  const Self& operator=(const Self& arc) {
    if (this == &arc) return (*this);

    // Free any existing data.
    if (m_extra_data != nullptr) delete m_extra_data;

    // Copy the arc's attributes.
    m_r = arc.m_r;
    m_s = arc.m_s;
    m_t = arc.m_t;
    m_u = arc.m_u;
    m_v = arc.m_v;
    m_w = arc.m_w;

    m_orient = arc.m_orient;
    m_info = arc.m_info;
    m_source = arc.m_source;
    m_target = arc.m_target;

    // Duplicate the extra data, if necessary.
    m_extra_data = (arc.m_extra_data != nullptr) ?
      new Extra_data(*(arc.m_extra_data)) : nullptr;

    return (*this);
  }
  //@}

  /// \name Get the arc properties.
  //@{

  /*! Check wheather the arc is valid.
   */
  bool is_valid() const { return test_flag(IS_VALID); }

  /*! Check whether the arc represents a full conic curve.
   */
  bool is_full_conic() const { return test_flag(IS_FULL_CONIC); }

  /*! Get the coefficients of the underlying conic.
   */
  const Integer& r() const { return (m_r); }
  const Integer& s() const { return (m_s); }
  const Integer& t() const { return (m_t); }
  const Integer& u() const { return (m_u); }
  const Integer& v() const { return (m_v); }
  const Integer& w() const { return (m_w); }

  /*! Check whether the arc is x-monotone.
   */
  bool is_x_monotone() const {
    // Check if the arc contains no vertical tangency points.
    Point_2 vtan_ps[2];
    return (vertical_tangency_points(vtan_ps) == 0);
  }

  /*! Check whether the arc is y-monotone.
   */
  bool is_y_monotone() const {
    // Check if the arc contains no horizontal tangency points.
    Point_2 htan_ps[2];
    return (horizontal_tangency_points(htan_ps) == 0);
  }

  /*! Obtain the arc's source.
   * \return The source point.
   * \pre The arc does not represent a full conic curve.
   */
  const Point_2& source() const {
    CGAL_precondition(! is_full_conic());
    return m_source;
  }

  /*! Obtain the arc's target.
   * \return The target point.
   * \pre The arc does not represent a full conic curve.
   */
  const Point_2& target() const {
    CGAL_precondition(! is_full_conic());
    return m_target;
  }

  /*! Obtain the orientation of the arc.
   * \return The orientation.
   */
  Orientation orientation() const { return m_orient; }

  /*! Obtain the extra data.
   */
  const Extra_data* extra_data() const { return m_extra_data; }

  /*! Obtain a bounding box for the conic arc.
   * \return The bounding box.
   */
  Bbox_2 bbox() const {
    CGAL_precondition(is_valid());

    double x_min(0), y_min(0), x_max(0), y_max(0);

    if (is_full_conic()) {
      // In case of a full conic (an ellipse or a circle), compute the
      // horizontal and vertical tangency points and use them to bound the arc.
      Point_2   tan_ps[2];
      CGAL_assertion_code(int n_tan_ps);

      CGAL_assertion_code(n_tan_ps = vertical_tangency_points(tan_ps));
      CGAL_assertion(n_tan_ps == 2);

      if (CGAL::to_double(tan_ps[0].x()) < CGAL::to_double(tan_ps[1].x())) {
        x_min = CGAL::to_double(tan_ps[0].x());
        x_max = CGAL::to_double(tan_ps[1].x());
      }
      else {
        x_min = CGAL::to_double(tan_ps[1].x());
        x_max = CGAL::to_double(tan_ps[0].x());
      }

      CGAL_assertion_code(n_tan_ps = horizontal_tangency_points(tan_ps));
      CGAL_assertion(n_tan_ps == 2);

      if (CGAL::to_double(tan_ps[0].y()) < CGAL::to_double(tan_ps[1].y())) {
        y_min = CGAL::to_double(tan_ps[0].y());
        y_max = CGAL::to_double(tan_ps[1].y());
      }
      else {
        y_min = CGAL::to_double(tan_ps[1].y());
        y_max = CGAL::to_double(tan_ps[0].y());
      }
    }
    else {
      // Use the source and target to initialize the exterme points.
      bool source_left =
        CGAL::to_double(m_source.x()) < CGAL::to_double(m_target.x());
      x_min = source_left ?
        CGAL::to_double(m_source.x()) : CGAL::to_double(m_target.x());
      x_max = source_left ?
        CGAL::to_double(m_target.x()) : CGAL::to_double(m_source.x());

      bool source_down =
        CGAL::to_double(m_source.y()) < CGAL::to_double(m_target.y());
      y_min = source_down ?
        CGAL::to_double(m_source.y()) : CGAL::to_double(m_target.y());
      y_max = source_down ?
        CGAL::to_double(m_target.y()) : CGAL::to_double(m_source.y());

      // Go over the vertical tangency points and try to update the x-points.
      Point_2 tan_ps[2];
      int n_tan_ps;

      n_tan_ps = vertical_tangency_points(tan_ps);
      for (int i = 0; i < n_tan_ps; ++i) {
        if (CGAL::to_double(tan_ps[i].x()) < x_min)
          x_min = CGAL::to_double(tan_ps[i].x());
        if (CGAL::to_double(tan_ps[i].x()) > x_max)
          x_max = CGAL::to_double(tan_ps[i].x());
      }

      // Go over the horizontal tangency points and try to update the y-points.
      n_tan_ps = horizontal_tangency_points(tan_ps);
      for (int i = 0; i < n_tan_ps; ++i) {
        if (CGAL::to_double(tan_ps[i].y()) < y_min)
          y_min = CGAL::to_double(tan_ps[i].y());
        if (CGAL::to_double(tan_ps[i].y()) > y_max)
          y_max = CGAL::to_double(tan_ps[i].y());
      }
    }

    // Return the resulting bounding box.
    return Bbox_2(x_min, y_min, x_max, y_max);
  }

  //@}

public:
  /// \name Flag manipulation functions.
  //@{
  template <typename T>
  constexpr size_t flag_mask(const T flag) const { return 0x1 << flag; }

  void reset_flags() { m_info = 0; }

  template <typename T>
  void set_flag(const T flag) { m_info |= flag_mask(flag); }

  template <typename T>
  void reset_flag(const T flag) { m_info &= ~flag_mask(flag); }

  template <typename T>
  void flip_flag(const T flag) { m_info ^= flag_mask(flag);}

  template <typename T>
  bool test_flag(const T flag) const
  { return (m_info & flag_mask(flag)); }
  //@}

public:
  /// \name Modifying functions (setters);
  // only friends have the priviledge to use.
  //@{

  /*! Set the source point of the conic arc.
   * \param ps The new source point.
   */
  void set_source(const Point_2& ps) { m_source = ps; }

  /*! Set the target point of the conic arc.
   * \param pt The new source point.
   */
  void set_target(const Point_2& pt) { m_target = pt; }

  /*! Set the coefficients.
   */
  void set_coefficients(Integer r, Integer s, Integer t,
                        Integer u, Integer v, Integer w) {

    m_r = r;
    m_s = s;
    m_t = t;
    m_u = u;
    m_v = v;
    m_w = w;
  }

  /*! Set the orientation.
   */
  void set_orientation(Orientation orient) { m_orient = orient; }

  /*! Set the endpoints.
   */
  void set_endpoints(const Conic_point_2& source, const Conic_point_2& target) {
    m_source = source;
    m_target = target;
  }

  /*! Set the extra data field.
   */
  void set_extra_data(Extra_data* extra_data) { m_extra_data = extra_data; }

  //@}

  /// \name Compute points on the arc.
  //@{

  /*! Calculate the vertical tangency points of the arc.
   * \param vpts The vertical tangency points.
   * \pre The vpts vector should be allocated at the size of 2.
   * \return The number of vertical tangency points.
   */
  int vertical_tangency_points(Point_2* vpts) const {
    // No vertical tangency points for line segments:
    if (m_orient == COLLINEAR) return 0;

    // Calculate the vertical tangency points of the supporting conic.
    Point_2 ps[2];
    auto n = conic_vertical_tangency_points(ps);

    // Return only the points that are contained in the arc interior.
    int m = 0;

    for (int i = 0; i < n; ++i) {
      if (is_full_conic() || is_strictly_between_endpoints(ps[i])) {
        vpts[m] = ps[i];
        ++m;
      }
    }

    // Return the number of vertical tangency points found.
    CGAL_assertion(m <= 2);
    return m;
  }

  /*! Calculate the horizontal tangency points of the arc.
   * \param hpts The horizontal tangency points.
   * \pre The hpts vector should be allocated at the size of 2.
   * \return The number of horizontal tangency points.
   */
  int horizontal_tangency_points(Point_2* hpts) const {
    // No horizontal tangency points for line segments:
    if (m_orient == COLLINEAR) return 0;

    // Calculate the horizontal tangency points of the conic.
    Point_2 ps[2];
    int n = conic_horizontal_tangency_points(ps);

    // Return only the points that are contained in the arc interior.
    int m = 0;

    for (int i = 0; i < n; ++i) {
      if (is_full_conic() || is_strictly_between_endpoints(ps[i])) {
        hpts[m] = ps[i];
        ++m;
      }
    }

    // Return the number of horizontal tangency points found.
    CGAL_assertion(m <= 2);
    return m;
  }

  /*! Find all points on the arc with a given x-coordinate.
   * \param p A placeholder for the x-coordinate.
   * \param ps The point on the arc at x(p).
   * \pre The vector ps should be allocated at the size of 2.
   * \return The number of points found.
   */
  int points_at_x(const Point_2& p, Point_2* ps) const {
    // Get the y coordinates of the points on the conic.
    Algebraic ys[2];
    int n = conic_get_y_coordinates(p.x(), ys);

    // Find all the points that are contained in the arc.
    int m = 0;

    for (int i = 0; i < n; ++i) {
      ps[m] = Point_2(p.x(), ys[i]);

      if (is_full_conic() || is_between_endpoints(ps[m])) ++m;
    }

    // Return the number of points on the arc.
    CGAL_assertion(m <= 2);
    return m;
  }

  /*! Find all points on the arc with a given y-coordinate.
   * \param p A placeholder for the y-coordinate.
   * \param ps The point on the arc at x(p).
   * \pre The vector ps should be allocated at the size of 2.
   * \return The number of points found.
   */
  int points_at_y(const Point_2& p, Point_2* ps) const {
    // Get the y coordinates of the points on the conic.
    Algebraic xs[2];
    int n = conic_get_x_coordinates(p.y(), xs);

    // Find all the points that are contained in the arc.
    int m = 0;

    for (int i = 0; i < n; ++i) {
      ps[m] = Point_2(xs[i], p.y());
      if (is_full_conic() || is_between_endpoints(ps[m])) ++m;
    }

    // Return the number of points on the arc.
    CGAL_assertion(m <= 2);
    return m;
  }
  //@}

private:
  /// \name Auxiliary construction functions.
  //@{

  /*! Set the properties of a conic arc (for the usage of the constructors).
   * \param rat_coeffs A vector of size 6, storing the rational coefficients
   *                   of x^2, y^2, xy, x, y and the free coefficient resp.
   */
  void set(const Rational* rat_coeffs) {
    set_flag(IS_VALID);

    // Convert the coefficients vector to an equivalent vector of integer
    // coefficients.
    Nt_traits nt_traits;
    Integer int_coeffs[6];

    nt_traits.convert_coefficients(rat_coeffs, rat_coeffs + 6, int_coeffs);

    // Check the orientation of conic curve, and negate the conic coefficients
    // if its given orientation.
    typename Rat_kernel::Conic_2 temp_conic(rat_coeffs[0], rat_coeffs[1],
                                            rat_coeffs[2], rat_coeffs[3],
                                            rat_coeffs[4], rat_coeffs[5]);

    if (m_orient == temp_conic.orientation()) {
      m_r = int_coeffs[0];
      m_s = int_coeffs[1];
      m_t = int_coeffs[2];
      m_u = int_coeffs[3];
      m_v = int_coeffs[4];
      m_w = int_coeffs[5];
    }
    else {
      m_r = -int_coeffs[0];
      m_s = -int_coeffs[1];
      m_t = -int_coeffs[2];
      m_u = -int_coeffs[3];
      m_v = -int_coeffs[4];
      m_w = -int_coeffs[5];
    }

    // Make sure both endpoint lie on the supporting conic.
    if (! is_on_supporting_conic(m_source) ||
        ! is_on_supporting_conic(m_target))
    {
      reset_flags();            // inavlid arc
      return;
    }

    m_extra_data = nullptr;

    // Check whether we have a degree 2 curve.
    if ((CGAL::sign(m_r) != ZERO) || (CGAL::sign(m_s) != ZERO) ||
        (CGAL::sign(m_t) != ZERO))
    {
      if (m_orient == COLLINEAR) {
        // We have a segment of a line pair with rational coefficients.
        // Compose the equation of the underlying line
        // (with algebraic coefficients).
        const Algebraic x1 = m_source.x(), y1 = m_source.y();
        const Algebraic x2 = m_target.x(), y2 = m_target.y();

        // The supporting line is A*x + B*y + C = 0, where:
        //
        //  A = y2 - y1,    B = x1 - x2,    C = x2*y1 - x1*y2
        //
        // We use the extra dat field to store the equation of this line.
        m_extra_data = new Extra_data;
        m_extra_data->a = y2 - y1;
        m_extra_data->b = x1 - x2;
        m_extra_data->c = x2*y1 - x1*y2;
        m_extra_data->side = ZERO;

        // Make sure the midpoint is on the line pair (thus making sure that
        // the two points are not taken from different lines).
        Alg_kernel ker;
        Point_2 p_mid = ker.construct_midpoint_2_object()(m_source, m_target);

        if (CGAL::sign((nt_traits.convert(m_r)*p_mid.x() +
                        nt_traits.convert(m_t)*p_mid.y() +
                        nt_traits.convert(m_u)) * p_mid.x() +
                       (nt_traits.convert(m_s)*p_mid.y() +
                        nt_traits.convert(m_v)) * p_mid.y() +
                       nt_traits.convert(m_w)) != ZERO)
        {
          reset_flags();        // inavlid arc
          return;
        }
      }
      else {
        // The sign of (4rs - t^2) detetmines the conic type:
        // - if it is possitive, the conic is an ellipse,
        // - if it is negative, the conic is a hyperbola,
        // - if it is zero, the conic is a parabola.
        CGAL::Sign sign_conic = CGAL::sign(4*m_r*m_s - m_t*m_t);

        // Build the extra hyperbolic data if necessary
        if (sign_conic == NEGATIVE) build_hyperbolic_arc_data();

        if (sign_conic != POSITIVE) {
          // In case of a non-degenerate parabola or a hyperbola, make sure
          // the arc is not infinite.
          Alg_kernel ker;
          Point_2 p_mid = ker.construct_midpoint_2_object()(m_source, m_target);
          Point_2 ps[2];

          bool finite_at_x = (points_at_x(p_mid, ps) > 0);
          bool finite_at_y = (points_at_y(p_mid, ps) > 0);

          if (! finite_at_x && ! finite_at_y) {
            reset_flags();              // inavlid arc
            return;
          }
        }
      }
    }


    set_flag(IS_VALID);         // mark that this arc valid
    reset_flag(IS_FULL_CONIC);  // mark that this arc is not a full conic
  }

  /*! Set the properties of a conic arc that is really a full curve
   * (that is, an ellipse).
   * \param rat_coeffs A vector of size 6, storing the rational coefficients
   *                   of x^2, y^2, xy, x, y and the free coefficient resp.
   * \param comp_orient Should we compute the orientation of the given curve.
   */
  void set_full(const Rational* rat_coeffs, const bool& comp_orient) {
    // Convert the coefficients vector to an equivalent vector of integer
    // coefficients.
    Nt_traits nt_traits;
    Integer int_coeffs[6];

    nt_traits.convert_coefficients(rat_coeffs, rat_coeffs + 6, int_coeffs);

    // Check the orientation of conic curve, and negate the conic coefficients
    // if its given orientation.
    typename Rat_kernel::Conic_2 temp_conic(rat_coeffs[0], rat_coeffs[1],
                                            rat_coeffs[2], rat_coeffs[3],
                                            rat_coeffs[4], rat_coeffs[5]);
    const Orientation temp_orient = temp_conic.orientation();

    if (comp_orient) m_orient = temp_orient;

    if (m_orient == temp_orient) {
      m_r = int_coeffs[0];
      m_s = int_coeffs[1];
      m_t = int_coeffs[2];
      m_u = int_coeffs[3];
      m_v = int_coeffs[4];
      m_w = int_coeffs[5];
    }
    else {
      m_r = -int_coeffs[0];
      m_s = -int_coeffs[1];
      m_t = -int_coeffs[2];
      m_u = -int_coeffs[3];
      m_v = -int_coeffs[4];
      m_w = -int_coeffs[5];
    }

    // Make sure the conic is a non-degenerate ellipse:
    // The coefficients should satisfy (4rs - t^2) > 0.
    const bool is_ellipse = (CGAL::sign(4*m_r*m_s - m_t*m_t) == POSITIVE);
    CGAL_assertion(is_ellipse);

    // We do not have to store any extra data with the arc.
    m_extra_data = nullptr;

    // Mark that this arc is a full conic curve.
    if (is_ellipse) {
      set_flag(IS_VALID);
      set_flag(IS_FULL_CONIC);
    }
    else reset_flags();            // inavlid arc
  }

  /*! Build the data for hyperbolic arc, contaning the characterization of the
   * hyperbolic branch the arc is placed on.
   */
  void build_hyperbolic_arc_data() {
    // Let phi be the rotation angle of the conic from its canonic form.
    // We can write:
    //
    //                          t
    //  sin(2*phi) = -----------------------
    //                sqrt((r - s)^2 + t^2)
    //
    //                        r - s
    //  cos(2*phi) = -----------------------
    //                sqrt((r - s)^2 + t^2)
    //
    Nt_traits nt_traits;
    const int or_fact = (m_orient == CLOCKWISE) ? -1 : 1;
    const Algebraic r = nt_traits.convert(or_fact * m_r);
    const Algebraic s = nt_traits.convert(or_fact * m_s);
    const Algebraic t = nt_traits.convert(or_fact * m_t);
    const Algebraic cos_2phi = (r - s) / nt_traits.sqrt((r-s)*(r-s) + t*t);
    const Algebraic zero = 0;
    const Algebraic one = 1;
    const Algebraic two = 2;
    Algebraic sin_phi;
    Algebraic cos_phi;

    // Calculate sin(phi) and cos(phi) according to the half-angle formulae:
    //
    //  sin(phi)^2 = 0.5 * (1 - cos(2*phi))
    //  cos(phi)^2 = 0.5 * (1 + cos(2*phi))
    Sign sign_t = CGAL::sign(t);

    if (sign_t == ZERO) {
      // sin(2*phi) == 0, so phi = 0 or phi = PI/2
      if (CGAL::sign(cos_2phi) == POSITIVE) {
        // phi = 0.
        sin_phi = zero;
        cos_phi = one;
      }
      else {
        // phi = PI/2.
        sin_phi = one;
        cos_phi = zero;
      }
    }
    else if (sign_t == POSITIVE) {
      // sin(2*phi) > 0 so 0 < phi < PI/2.
      sin_phi = nt_traits.sqrt((one + cos_2phi) / two);
      cos_phi = nt_traits.sqrt((one - cos_2phi) / two);
    }
    else {
      // sin(2*phi) < 0 so PI/2 < phi < PI.
      sin_phi = nt_traits.sqrt((one + cos_2phi) / two);
      cos_phi = -nt_traits.sqrt((one - cos_2phi) / two);
    }

    // Calculate the center (x0, y0) of the conic, given by the formulae:
    //
    //        t*v - 2*s*u                t*u - 2*r*v
    //  x0 = -------------   ,     y0 = -------------
    //        4*r*s - t^2                4*r*s - t^2
    //
    // The denominator (4*r*s - t^2) must be negative for hyperbolas.
    const Algebraic u = nt_traits.convert(or_fact * m_u);
    const Algebraic v = nt_traits.convert(or_fact * m_v);
    const Algebraic det = 4*r*s - t*t;
    Algebraic x0, y0;

    CGAL_assertion(CGAL::sign(det) == NEGATIVE);

    x0 = (t*v - two*s*u) / det;
    y0 = (t*u - two*r*v) / det;

    // The axis separating the two branches of the hyperbola is now given by:
    //
    //  cos(phi)*x + sin(phi)*y - (cos(phi)*x0 + sin(phi)*y0) = 0
    //
    // We store the equation of this line in the extra data structure and also
    // the sign (side of half-plane) our arc occupies with respect to the line.
    m_extra_data = new Extra_data;

    m_extra_data->a = cos_phi;
    m_extra_data->b = sin_phi;
    m_extra_data->c = - (cos_phi*x0 + sin_phi*y0);

    // Make sure that the two endpoints are located on the same branch
    // of the hyperbola.
    m_extra_data->side = sign_of_extra_data(m_source.x(), m_source.y());

    CGAL_assertion(m_extra_data->side != ZERO);
    CGAL_assertion(m_extra_data->side ==
                   sign_of_extra_data(m_target.x(), m_target.y()));
  }
  //@}

public:
  /// \name Auxiliary functions.
  //@{

  /*! Evaluate the sign of (a*x + b*y + c) stored with the extra data field
   * at a given point.
   * \param px The x-coordinate of query point.
   * \param py The y-coordinate of query point.
   * \return The sign of (a*x + b*y + c).
   */
  Sign sign_of_extra_data(const Algebraic& px, const Algebraic& py) const {
    CGAL_assertion(m_extra_data != nullptr);
    if (m_extra_data == nullptr) return ZERO;
    Algebraic val = m_extra_data->a*px + m_extra_data->b*py + m_extra_data->c;
    return CGAL::sign(val);
  }

protected:
  /*! Check whether the given point lies on the supporting conic of the arc.
   * \param p The query point.
   * \return true if p lies on the supporting conic; (false) otherwise.
   */
  bool is_on_supporting_conic(const Point_2& p) const {
    // Check whether p satisfies the conic equation.
    // The point must satisfy: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0.
    Nt_traits nt_traits;
    const Algebraic val = (nt_traits.convert(m_r)*p.x() +
                           nt_traits.convert(m_t)*p.y() +
                           nt_traits.convert(m_u)) * p.x() +
      (nt_traits.convert(m_s) * p.y() + nt_traits.convert(m_v)) * p.y() +
      nt_traits.convert(m_w);

    return (CGAL::sign(val) == ZERO);
  }

  /*! Check whether the given point is between the source and the target.
   * The point is assumed to be on the conic's boundary.
   * \param p The query point.
   * \return true if the point is between the two endpoints,
   *         (false) if it is not.
   */
  bool is_between_endpoints(const Point_2& p) const {
    CGAL_precondition(! is_full_conic());

    // Check if p is one of the endpoints.
    Alg_kernel ker;

    if (ker.equal_2_object()(p, m_source) || ker.equal_2_object()(p, m_target))
      return true;
    else return (is_strictly_between_endpoints(p));
  }

  /*! Check whether the given point is strictly between the source and the
   * target (but not any of them).
   * The point is assumed to be on the conic's boundary.
   * \param p The query point.
   * \return true if the point is strictly between the two endpoints,
   *         (false) if it is not.
   */
  bool is_strictly_between_endpoints(const Point_2& p) const {
    // In case this is a full conic, any point on its boundary is between
    // its end points.
    if (is_full_conic()) return true;

    // Check if we have extra data available.
    if (m_extra_data != nullptr) {
      if (m_extra_data->side != ZERO) {
        // In case of a hyperbolic arc, make sure the point is located on the
        // same branch as the arc.
        if (sign_of_extra_data(p.x(), p.y()) != m_extra_data->side)
          return false;
      }
      else {
        // In case we have a segment of a line pair, make sure that p really
        // satisfies the equation of the line.
        if (sign_of_extra_data(p.x(), p.y()) != ZERO) return false;
      }
    }

    // Act according to the conic degree.
    Alg_kernel ker;

    if (m_orient == COLLINEAR) {
      Comparison_result res1;
      Comparison_result res2;

      if (ker.compare_x_2_object()(m_source, m_target) == EQUAL) {
        // In case of a vertical segment - just check whether the y coordinate
        // of p is between those of the source's and of the target's.
        res1 = ker.compare_y_2_object()(p, m_source);
        res2 = ker.compare_y_2_object()(p, m_target);
      }
      else {
        // Otherwise, since the segment is x-monotone, just check whether the
        // x coordinate of p is between those of the source's and of the
        // target's.
        res1 = ker.compare_x_2_object()(p, m_source);
        res2 = ker.compare_x_2_object()(p, m_target);
      }

      // If p is not in the (open) x-range (or y-range) of the segment, it
      // cannot be contained in the segment.
      if ((res1 == EQUAL) || (res2 == EQUAL) || (res1 == res2)) return false;

      // Perform an orientation test: This is crucial for segment of line
      // pairs, as we want to make sure that p lies on the same line as the
      // source and the target.
      return (ker.orientation_2_object()(m_source, p, m_target) == COLLINEAR);
    }
    else {
      // In case of a conic of degree 2, make a decision based on the conic's
      // orientation and whether (source,p,target) is a right or a left turn.
      if (m_orient == COUNTERCLOCKWISE)
        return (ker.orientation_2_object()(m_source, p, m_target) == LEFT_TURN);
      else
        return (ker.orientation_2_object()(m_source, p, m_target) == RIGHT_TURN);
    }
  }

  /*! Find the vertical tangency points of the undelying conic.
   * \param ps The output points of vertical tangency.
   *           This area must be allocated at the size of 2.
   * \return The number of vertical tangency points.
   */
  int conic_vertical_tangency_points(Point_2* ps) const {
    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    if (CGAL::sign(m_s) == ZERO) return 0;

    // We are interested in the x coordinates where the quadratic equation:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    // has a single solution (obviously if s = 0, there are no such points).
    // We therefore demand that the discriminant of this equation is zero:
    //  (t*x + v)^2 - 4*s*(r*x^2 + u*x + w) = 0
    const Integer two = 2;
    const Integer four = 4;
    Algebraic xs[2];
    Algebraic* xs_end;
    int n_xs;
    Nt_traits nt_traits;

    xs_end = nt_traits.solve_quadratic_equation(m_t*m_t - four*m_r*m_s,
                                                two*m_t*m_v - four*m_s*m_u,
                                                m_v*m_v - four*m_s*m_w,
                                                xs);
    n_xs = static_cast<int>(xs_end - xs);

    // Find the y-coordinates of the vertical tangency points.
    Algebraic ys[2];
    Algebraic* ys_end;
    int n_ys;

    if (CGAL::sign(m_t) == ZERO) {
      // The two vertical tangency points have the same y coordinate:
      ys[0] = nt_traits.convert(-m_v) / nt_traits.convert(two*m_s);
      n_ys = 1;
    }
    else {
      ys_end =
        nt_traits.solve_quadratic_equation(four*m_r*m_s*m_s - m_s*m_t*m_t,
                                           four*m_r*m_s*m_v - two*m_s*m_t*m_u,
                                           m_r*m_v*m_v - m_t*m_u*m_v + m_t*m_t*m_w,
                                           ys);
      n_ys = static_cast<int>(ys_end - ys);
    }

    // Pair the x and y coordinates and obtain the vertical tangency points.
    int n(0);

    for (int i = 0; i < n_xs; ++i) {
      if (n_ys == 1) {
        ps[n] = Point_2(xs[i], ys[0]);
        ++n;
      }
      else {
        for (int j = 0; j < n_ys; j++) {
          if (CGAL::compare(nt_traits.convert(two*m_s) * ys[j],
                            -(nt_traits.convert(m_t) * xs[i] +
                              nt_traits.convert(m_v))) == EQUAL)
          {
            ps[n] = Point_2(xs[i], ys[j]);
            ++n;
            break;
          }
        }
      }
    }

    CGAL_assertion(n <= 2);
    return n;
  }

  /*! Find the horizontal tangency points of the undelying conic.
   * \param ps The output points of horizontal tangency.
   *           This area must be allocated at the size of 2.
   * \return The number of horizontal tangency points.
   */
  int conic_horizontal_tangency_points(Point_2* ps) const {
    const Integer zero = 0;

    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    if (CGAL::sign(m_r) == ZERO) return 0;

    // We are interested in the y coordinates were the quadratic equation:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    // has a single solution (obviously if r = 0, there are no such points).
    // We therefore demand that the discriminant of this equation is zero:
    //  (t*y + u)^2 - 4*r*(s*y^2 + v*y + w) = 0
    const Integer two = 2;
    const Integer four = 4;
    int n;
    Algebraic ys[2];
    Algebraic *ys_end;
    Nt_traits nt_traits;

    ys_end = nt_traits.solve_quadratic_equation(m_t*m_t - four*m_r*m_s,
                                                two*m_t*m_u - four*m_r*m_v,
                                                m_u*m_u - four*m_r*m_w,
                                                ys);
    n = static_cast<int>(ys_end - ys);

    // Compute the x coordinates and construct the horizontal tangency points.
    for (int i = 0; i < n; ++i) {
      // Having computed y, x is the single solution to the quadratic equation
      // above, and since its discriminant is 0, x is simply given by:
      Algebraic x = -(nt_traits.convert(m_t)*ys[i] + nt_traits.convert(m_u)) /
        nt_traits.convert(two*m_r);
      ps[i] = Point_2(x, ys[i]);
    }

    CGAL_assertion(n <= 2);
    return n;
  }

  /*! Find the y coordinates of the underlying conic at a given x coordinate.
   * \param x The x coordinate.
   * \param ys The output y coordinates.
   * \pre The vector ys must be allocated at the size of 2.
   * \return The number of y coordinates computed (either 0, 1 or 2).
   */
  int conic_get_y_coordinates(const Algebraic& x, Algebraic* ys) const {
    // Solve the quadratic equation for a given x and find the y values:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    Nt_traits nt_traits;
    Algebraic A = nt_traits.convert(m_s);
    Algebraic B = nt_traits.convert(m_t)*x + nt_traits.convert(m_v);
    Algebraic C = (nt_traits.convert(m_r)*x + nt_traits.convert(m_u))*x +
      nt_traits.convert(m_w);

    return (solve_quadratic_equation(A, B, C, ys[0], ys[1]));
  }

  /*! Find the x coordinates of the underlying conic at a given y coordinate.
   * \param y The y coordinate.
   * \param xs The output x coordinates.
   * \pre The vector xs must be allocated at the size of 2.
   * \return The number of x coordinates computed (either 0, 1 or 2).
   */
  int conic_get_x_coordinates(const Algebraic& y, Algebraic* xs) const {
    // Solve the quadratic equation for a given y and find the x values:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    Nt_traits nt_traits;
    Algebraic A = nt_traits.convert(m_r);
    Algebraic B = nt_traits.convert(m_t)*y + nt_traits.convert(m_u);
    Algebraic C = (nt_traits.convert(m_s)*y + nt_traits.convert(m_v))*y +
      nt_traits.convert(m_w);

    return (solve_quadratic_equation(A, B, C, xs[0], xs[1]));
  }

  /*! Solve the given quadratic equation: Ax^2 + B*x + C = 0.
   * \param x_minus The root obtained from taking -sqrt(discriminant).
   * \param x_plus The root obtained from taking -sqrt(discriminant).
   * \return The number of disticnt solutions to the equation.
   */
  int solve_quadratic_equation(const Algebraic& A,
                               const Algebraic& B,
                               const Algebraic& C,
                               Algebraic& x_minus, Algebraic& x_plus) const {
    // Check if we actually have a linear equation.
    if (CGAL::sign(A) == ZERO) {
      if (CGAL::sign(B) == ZERO) return 0;
      x_minus = x_plus = -C / B;
      return 1;
    }

    // Compute the discriminant and act according to its sign.
    const Algebraic disc = B*B - 4*A*C;
    Sign sign_disc = CGAL::sign(disc);

    // Check whether there are no real-valued solutions:
    if (sign_disc == NEGATIVE) return 0;
    else if (sign_disc == ZERO) {
      // One distinct solution:
      x_minus = x_plus = -B / (2*A);
      return 1;
    }

    // Compute the two distinct solutions:
    Algebraic _2A = 2*A;
    Nt_traits nt_traits;
    Algebraic sqrt_disc = nt_traits.sqrt (disc);

    x_minus = -(B + sqrt_disc) / _2A;
    x_plus = (sqrt_disc - B) / _2A;
    return 2;
  }
  //@}

};

/*! Exporter for conic arcs.
 */
template <typename Rat_kernel, typename Alg_kernel, typename Nt_traits>
std::ostream&
operator<< (std::ostream& os,
            const _Conic_arc_2<Rat_kernel, Alg_kernel, Nt_traits>& arc)
{
  os << "{" << CGAL::to_double(arc.r()) << "*x^2 + "
     << CGAL::to_double(arc.s()) << "*y^2 + "
     << CGAL::to_double(arc.t()) << "*xy + "
     << CGAL::to_double(arc.u()) << "*x + "
     << CGAL::to_double(arc.v()) << "*y + "
     << CGAL::to_double(arc.w()) << "}";

  if (arc.is_full_conic()) os << " - Full curve";
  else {
    os << " : (" << CGAL::to_double(arc.source().x()) << ","
       << CGAL::to_double(arc.source().y()) << ") ";

    if (arc.orientation() == CLOCKWISE)
      os << "--cw-->";
    else if (arc.orientation() == COUNTERCLOCKWISE)
      os << "--ccw-->";
    else
      os << "--l-->";

    os << " (" << CGAL::to_double(arc.target().x()) << ","
       << CGAL::to_double(arc.target().y()) << ")";
  }

  return (os);
}

} //namespace CGAL

#endif
