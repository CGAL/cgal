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

  /*! Determine wheather the arc is valid.
   */
  bool is_valid() const { return test_flag(IS_VALID); }

  /*! Determine whether the arc represents a full conic curve.
   */
  bool is_full_conic() const { return test_flag(IS_FULL_CONIC); }

  /*! Obtain the coefficients of the underlying conic.
   */
  const Integer& r() const { return (m_r); }
  const Integer& s() const { return (m_s); }
  const Integer& t() const { return (m_t); }
  const Integer& u() const { return (m_u); }
  const Integer& v() const { return (m_v); }
  const Integer& w() const { return (m_w); }

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

protected:
  template <typename, typename, typename> friend class Arr_conic_traits_2;

  /// \name Flag manipulation functions.
  //@{
  template <typename T>
  static constexpr size_t flag_mask(const T flag) { return 0x1 << flag; }

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
        for (int j = 0; j < n_ys; ++j) {
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
    const Integer zero(0);

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
    Algebraic* ys_end;
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
