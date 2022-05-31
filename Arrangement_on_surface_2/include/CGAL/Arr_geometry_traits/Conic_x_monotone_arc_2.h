// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_CONIC_X_MONOTONE_ARC_2_H
#define CGAL_CONIC_X_MONOTONE_ARC_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Header file for the _Conic_x_monotone_arc_2<Conic_arc_2> class.
 */

#include <map>
#include <ostream>

#include <boost/variant.hpp>

#include <CGAL/Arr_geometry_traits/Conic_intersections_2.h>

namespace CGAL {

/*! Representation of an x-monotone conic arc.
 * The class is templated by a representation of a general bounded conic arc.
 */

template <typename ConicArc>
class _Conic_x_monotone_arc_2 : public ConicArc {
public:
  typedef ConicArc                                Conic_arc_2;
  typedef _Conic_x_monotone_arc_2<Conic_arc_2>    Self;

  typedef typename Conic_arc_2::Alg_kernel        Alg_kernel;
  typedef typename Conic_arc_2::Algebraic         Algebraic;

  typedef typename Conic_arc_2::Point_2           Point_2;
  typedef typename Conic_arc_2::Conic_point_2     Conic_point_2;

  // Type definition for the intersection points mapping.
  typedef typename Conic_point_2::Conic_id        Conic_id;

  using Conic_arc_2::sign_of_extra_data;

  using Conic_arc_2::IS_VALID;
  using Conic_arc_2::IS_FULL_CONIC;

  using Conic_arc_2::flag_mask;
  using Conic_arc_2::reset_flags;
  using Conic_arc_2::set_flag;
  using Conic_arc_2::reset_flag;
  using Conic_arc_2::flip_flag;
  using Conic_arc_2::test_flag;

private:
  typedef Conic_arc_2                             Base;

  typedef typename Conic_arc_2::Integer           Integer;
  typedef typename Conic_arc_2::Nt_traits         Nt_traits;
  typedef typename Conic_arc_2::Rat_kernel        Rat_kernel;

  Algebraic m_alg_r;    // The coefficients of the supporting conic curve:
  Algebraic m_alg_s;    //
  Algebraic m_alg_t;    //   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0 ,
  Algebraic m_alg_u;    //
  Algebraic m_alg_v;    // converted to algebraic numbers.
  Algebraic m_alg_w;    //

  Conic_id m_id;        // The ID number of the supporting conic curve.

public:
  // Bit masks for the m_info field (the two least significant bits are already
  // used by the base class).
  enum {
    IS_VERTICAL_SEGMENT = Conic_arc_2::LAST_INFO,
    IS_DIRECTED_RIGHT,
    DEGREE_1,
    DEGREE_2,
    PLUS_SQRT_DISC_ROOT,
    FACING_UP,
    FACING_DOWN,
    IS_SPECIAL_SEGMENT,

    DEGREE_MASK = (0x1 << DEGREE_1) | (0x1 << DEGREE_2),
    FACING_MASK = (0x1 << FACING_UP) | (0x1 << FACING_DOWN)
  };

  /// \name Public constrcutors, assignment operators, and destructors.
  //@{

  /*! Default constructor.
   */
  _Conic_x_monotone_arc_2() :
    Base(),
    m_id()
  {}

  /*! Copy constructor.
   * \param arc The copied arc.
   */
  _Conic_x_monotone_arc_2(const Self& arc) :
    Base(arc),
    m_alg_r(arc.m_alg_r),
    m_alg_s(arc.m_alg_s),
    m_alg_t(arc.m_alg_t),
    m_alg_u(arc.m_alg_u),
    m_alg_v(arc.m_alg_v),
    m_alg_w(arc.m_alg_w),
    m_id(arc.m_id)
  {}

  /*! Assignment operator.
   * \param arc The copied arc.
   */
  const Self& operator=(const Self& arc) {
    CGAL_precondition (arc.is_valid());

    if (this == &arc) return (*this);

    // Copy the base arc.
    Base::operator= (arc);

    // Set the rest of the properties.
    m_alg_r = arc.m_alg_r;
    m_alg_s = arc.m_alg_s;
    m_alg_t = arc.m_alg_t;
    m_alg_u = arc.m_alg_u;
    m_alg_v = arc.m_alg_v;
    m_alg_w = arc.m_alg_w;

    m_id = arc.m_id;

    return (*this);
  }
  //@}

private:
  template <typename, typename, typename> friend class Arr_conic_traits_2;

  /// \name private constrcutors to be used only by the traits class template.
  //@{

  /*! Construct an x-monotone arc from a conic arc.
   * \param arc The given (base) arc.
   * \pre The given arc is x-monotone.
   */
  _Conic_x_monotone_arc_2(const Base& arc) : Base(arc), m_id() {}

  /*! Construct an x-monotone arc from a conic arc.
   * \param arc The given (base) arc.
   * \param id The ID of the base arc.
   */
  _Conic_x_monotone_arc_2(const Base& arc, const Conic_id& id) :
    Base(arc),
    m_id(id)
  {}

  /*! Construct an x-monotone sub-arc from a conic arc.
   * \param arc The given (base) arc.
   * \param source The source point.
   * \param target The target point.
   * \param id The ID of the base arc.
   */
  _Conic_x_monotone_arc_2(const Base& arc,
                          const Point_2& source, const Point_2& target,
                          const Conic_id& id) :
    Base(arc),
    m_id(id)
  {}

  /*! Construct a special segment connecting to given endpoints (for the usage
   * of the landmarks point-location strategy).
   * \param source The source point.
   * \param target The target point.
   */
  _Conic_x_monotone_arc_2(const Point_2& source, const Point_2& target) :
    Base(source, target)
  {}

  /*! Construct a special segment of a given line connecting to given
   * endpoints.
   * \param a, b, c The coefficients of the supporting line (ax + by + c = 0).
   * \param source The source point.
   * \param target The target point.
   */
  _Conic_x_monotone_arc_2(const Algebraic& a,
                          const Algebraic& b,
                          const Algebraic& c,
                          const Point_2& source, const Point_2& target) :
    Base()
  {}

  /*! Obtain the non const reference to the curve source point.
   */
  Conic_point_2& source() { return (this->m_source); }

  /*! Obtain the non const reference to the curve source point.
   */
  Conic_point_2& target() { return this->m_target; }
  //@}

public:
  /// \name Accessing the arc properties.
  //@{

  /*! Obtain the facing mask.
   */
  size_t facing_mask() const { return this->m_info & FACING_MASK; }

  /*! Obtain the degree mask.
   */
  size_t degree_mask() const { return this->m_info & DEGREE_MASK; }

  /*! Obtain the coefficients of the underlying conic.
   */
  const Integer& r() const { return ( this->m_r); }
  const Integer& s() const { return ( this->m_s); }
  const Integer& t() const { return ( this->m_t); }
  const Integer& u() const { return ( this->m_u); }
  const Integer& v() const { return ( this->m_v); }
  const Integer& w() const { return ( this->m_w); }

  /*! Obtain the arc's source.
   * \return The source point.
   */
  const Conic_point_2& source() const { return (this->m_source); }

  /*! Obtain the arc's target.
   * \return The target point.
   */
  const Conic_point_2& target() const { return this->m_target; }

  /*! Obtain the orientation of the arc.
   * \return The orientation.
   */
  Orientation orientation() const { return this->m_orient; }

  /*! Obtain the left endpoint of the arc.
   */
  const Conic_point_2& left() const {
    if (test_flag(IS_DIRECTED_RIGHT)) return this->m_source;
    else return this->m_target;
  }

  /*! Obtain the right endpoint of the arc.
   */
  const Conic_point_2& right() const {
    if (test_flag(IS_DIRECTED_RIGHT)) return this->m_target;
    else return this->m_source;
  }

  /*! Determine whether the conic arc is directed iexicographically right.
   */
  bool is_directed_right() const { return test_flag(IS_DIRECTED_RIGHT); }

  /*! Determine whether the conic arc is a vertical segment.
   */
  bool is_vertical() const { return test_flag(IS_VERTICAL_SEGMENT); }

  /*! Determine whether the conic arc is a facing up.
   */
  bool is_upper() const { return test_flag(FACING_UP); }

  /*! Determine whether the conic arc is a facing down.
   */
  bool is_lower() const { return test_flag(FACING_DOWN); }

  /*! Check whether the arc is a special segment connecting two algebraic
   * endpoints (and has no undelying integer conic coefficients).
   */
  bool is_special_segment() const { return test_flag(IS_SPECIAL_SEGMENT); }

  /*! Obtain the mask of the DEGREE_1 flag.
   */
  static constexpr size_t degree_1_mask() { return flag_mask(DEGREE_1); }

  /*! Obtain the mask of the DEGREE_1 flag.
   */
  static constexpr size_t degree_2_mask() { return flag_mask(DEGREE_2); }

  /*! Obtain the algebraic coefficients.
   */
  Algebraic alg_r() const { return m_alg_r; }
  Algebraic alg_s() const { return m_alg_s; }
  Algebraic alg_t() const { return m_alg_t; }
  Algebraic alg_u() const { return m_alg_u; }
  Algebraic alg_v() const { return m_alg_v; }
  Algebraic alg_w() const { return m_alg_w; }

  /*! Obtain the conic id.
   */
  Conic_id id() const { return m_id; }

  /*! Obtain a bounding box for the conic arc.
   * \return The bounding box.
   */
  Bbox_2 bbox() const { return Base::bbox(); }
  //@}

  // Setters
  //@{

  /*! Set the algebraic coefficients.
   */
  void set_alg_coefficients(const Algebraic& alg_r, const Algebraic& alg_s,
                            const Algebraic& alg_t, const Algebraic& alg_u,
                            const Algebraic& alg_v, const Algebraic& alg_w)
  {
    m_alg_r = alg_r;
    m_alg_s = alg_s;
    m_alg_t = alg_t;
    m_alg_u = alg_u;
    m_alg_v = alg_v;
    m_alg_w = alg_w;
  }

  /*! Add a generating conic ID.
   */
  void set_generating_conic(const Conic_id& id) {
    this->m_source.set_generating_conic(id);
    this->m_target.set_generating_conic(id);
  }
  //@}

  /// \name Constructing points on the arc.
  //@{

  /*! Obtain a polyline approximating the conic arc.
   * \param n The maximal number of sample points.
   * \param oi An output iterator, whose value-type is pair<double,double>
   *           (representing an approximated point).
   *           In case the arc is a line segment, there are 2 output points,
   *           otherwise the arc is approximated by the polyline defined by
   *           (p_0, p_1, ..., p_n), where p_0 and p_n are the left and right
   *           endpoints of the arc, respectively.
   */
  template <typename OutputIterator>
  OutputIterator polyline_approximation(size_t n, OutputIterator oi) const {
    CGAL_precondition (n != 0);

    const double x_left = CGAL::to_double (left().x());
    const double y_left = CGAL::to_double (left().y());
    const double x_right = CGAL::to_double (right().x());
    const double y_right = CGAL::to_double (right().y());

    if (this->m_orient == COLLINEAR) {
      // In case of a line segment, return the two endpoints.
      *oi++ = std::pair<double, double>(x_left, y_left);
      *oi++ = std::pair<double, double>(x_right, y_right);
      return oi;
    }

    // Otherwise, sample (n - 1) equally-spaced points in between.
    const double app_r = CGAL::to_double(this->m_r);
    const double app_s = CGAL::to_double(this->m_s);
    const double app_t = CGAL::to_double(this->m_t);
    const double app_u = CGAL::to_double(this->m_u);
    const double app_v = CGAL::to_double(this->m_v);
    const double app_w = CGAL::to_double(this->m_w);
    const double x_jump = (x_right - x_left) / n;
    double x, y;
    const bool A_is_zero = (CGAL::sign(this->m_s) == ZERO);
    double A = app_s, B, C;
    double disc;
    size_t i;

    *oi++ = std::pair<double, double>(x_left, y_left);      // the left point
    for (i = 1; i < n; ++i) {
      x = x_left + x_jump*i;

      // Solve the quadratic equation: A*x^2 + B*x + C = 0:
      B = app_t*x + app_v;
      C = (app_r*x + app_u)*x + app_w;

      if (A_is_zero) y = -C / B;
      else {
        disc = B*B - 4*A*C;
        if (disc < 0) disc = 0;

        // We take either the root involving -sqrt(disc) or +sqrt(disc)
        // based on the information flags.
        y = (test_flag(PLUS_SQRT_DISC_ROOT)) ?
          (std::sqrt(disc) - B) / (2*A) : -(B + std::sqrt (disc)) / (2*A);
      }
      *oi++ = std::pair<double, double>(x, y);
    }
    *oi++ = std::pair<double, double>(x_right, y_right);    // the right point

    return oi;
  }

  //@}

  /// \name Constructing x-monotone arcs.
  //@{

  /*! Flip the arc.
   * \return An arc with swapped source and target and a reverse orienation.
   */
  Self flip() const {
    // Make a copy of the current arc.
    Self arc = *this;

    // Reverse the orientation.
    if (this->m_orient == CLOCKWISE) arc.m_orient = COUNTERCLOCKWISE;
    else if (this->m_orient == COUNTERCLOCKWISE) arc.m_orient = CLOCKWISE;

    // Swap the source and the target.
    arc.m_source = this->m_target;
    arc.m_target = this->m_source;

    // Change the direction bit among the information flags.
    arc.flip_flag(IS_DIRECTED_RIGHT);

    return arc;
  }
  //@}

private:
  /// \name Auxiliary (private) functions.
  //@{

  /*! Check whether the given point lies on the supporting conic of the arc.
   * \param px The x-coordinate of query point.
   * \param py The y-coordinate of query point.
   * \return (true) if p lies on the supporting conic; (false) otherwise.
   */
  bool is_on_supporting_conic(const Algebraic& px, const Algebraic& py) const {
    CGAL::Sign _sign;

    if (! is_special_segment()) {
      // Check whether p satisfies the conic equation.
      // The point must satisfy: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0.
      _sign = CGAL::sign((m_alg_r*px + m_alg_t*py + m_alg_u) * px +
                         (m_alg_s*py + m_alg_v) * py + m_alg_w);
    }
    else {
      // Check whether p satisfies the equation of the line stored with the
      // extra data.
      _sign = sign_of_extra_data(px, py);
    }

    return (_sign == ZERO);
  }

public:
  /*! Obtain the i'th order derivative by x of the conic at the point p=(x,y).
   * \param p The point where we derive.
   * \param i The order of the derivatives (either 1, 2 or 3).
   * \param slope_numer The numerator of the slope.
   * \param slope_denom The denominator of the slope.
   * \todo Allow higher order derivatives.
   */
  void derive_by_x_at(const Point_2& p, const unsigned int& i,
                      Algebraic& slope_numer, Algebraic& slope_denom) const {
    if (is_special_segment()) {
      // Special treatment for special segments, given by (a*x + b*y + c = 0),
      // so their first-order derivative by x is simply -a/b. The higher-order
      // derivatives are all 0.
      if (i == 1) {
        if (CGAL::sign (this->m_extra_data->b) != NEGATIVE) {
          slope_numer = - this->m_extra_data->a;
          slope_denom = this->m_extra_data->b;
        }
        else {
          slope_numer = this->m_extra_data->a;
          slope_denom = - this->m_extra_data->b;
        }
      }
      else {
        slope_numer = 0;
        slope_denom = 1;
      }

      return;
    }

    // The derivative by x of the conic
    //   C: {r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0}
    // at the point p=(x,y) is given by:
    //
    //           2r*x + t*y + u       alpha
    //   y' = - ---------------- = - -------
    //           2s*y + t*x + v       beta
    //
    const Algebraic two = 2;
    const Algebraic sl_numer = two*m_alg_r*p.x() + m_alg_t*p.y() + m_alg_u;
    const Algebraic sl_denom = two*m_alg_s*p.y() + m_alg_t*p.x() + m_alg_v;

    if (i == 1) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE) {
        slope_numer = -sl_numer;
        slope_denom = sl_denom;
      }
      else {
        slope_numer = sl_numer;
        slope_denom = -sl_denom;
      }

      return;
    }

    // The second-order derivative is given by:
    //
    //             s*alpha^2 - t*alpha*beta + r*beta^2     gamma
    //   y'' = -2 ------------------------------------- = -------
    //                           beta^3                    delta
    //
    const Algebraic sl2_numer = m_alg_s * sl_numer*sl_numer -
      m_alg_t * sl_numer*sl_denom + m_alg_r * sl_denom*sl_denom;
    const Algebraic sl2_denom = sl_denom*sl_denom*sl_denom;

    if (i == 2) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE) {
        slope_numer = -two *sl2_numer;
        slope_denom = sl2_denom;
      }
      else {
        slope_numer = two *sl2_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // The third-order derivative is given by:
    //
    //              (2s*alpha - t*beta) * gamma
    //   y''' = -6 ------------------------------
    //                    beta^2 * delta
    //
    const Algebraic sl3_numer =
      (two * m_alg_s * sl_numer - m_alg_t * sl_denom) * sl2_numer;
    const Algebraic sl3_denom = sl_denom*sl_denom * sl2_denom;

    if (i == 3) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE) {
        slope_numer = -6 * sl3_numer;
        slope_denom = sl3_denom;
      }
      else {
        slope_numer = 6 * sl3_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // \todo Handle higher-order derivatives as well.
    CGAL_error();
  }

  /*! Obtain the i'th order derivative by y of the conic at the point p=(x,y).
   * \param p The point where we derive.
   * \param i The order of the derivatives (either 1, 2 or 3).
   * \param slope_numer The numerator of the slope.
   * \param slope_denom The denominator of the slope.
   * \todo Allow higher order derivatives.
   */
  void derive_by_y_at(const Point_2& p, const int& i,
                      Algebraic& slope_numer, Algebraic& slope_denom) const {
    if (is_special_segment()) {
      // Special treatment for special segments, given by (a*x + b*y + c = 0),
      // so their first-order derivative by x is simply -b/a. The higher-order
      // derivatives are all 0.
      if (i == 1) {
        if (CGAL::sign(this->m_extra_data->a) != NEGATIVE) {
          slope_numer = - this->m_extra_data->b;
          slope_denom = this->m_extra_data->a;
        }
        else {
          slope_numer = this->m_extra_data->b;
          slope_denom = - this->m_extra_data->a;
        }
      }
      else {
        slope_numer = 0;
        slope_denom = 1;
      }

      return;
    }

    // The derivative by y of the conic
    //   C: {r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0}
    // at the point p=(x,y) is given by:
    //
    //           2s*y + t*x + v     alpha
    //   x' = - ---------------- = -------
    //           2r*x + t*y + u      beta
    //
    const Algebraic two = 2;
    const Algebraic sl_numer = two*m_alg_s*p.y() + m_alg_t*p.x() + m_alg_v;
    const Algebraic sl_denom = two*m_alg_r*p.x() + m_alg_t*p.y() + m_alg_u;

    if (i == 1) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE) {
        slope_numer = -sl_numer;
        slope_denom = sl_denom;
      }
      else {
        slope_numer = sl_numer;
        slope_denom = -sl_denom;
      }

      return;
    }

    // The second-order derivative is given by:
    //
    //             r*alpha^2 - t*alpha*beta + s*beta^2
    //   x'' = -2 -------------------------------------
    //                           beta^3
    //
    const Algebraic sl2_numer = m_alg_r * sl_numer*sl_numer -
      m_alg_t * sl_numer*sl_denom + m_alg_s * sl_denom*sl_denom;
    const Algebraic sl2_denom = sl_denom*sl_denom*sl_denom;

    if (i == 2) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE) {
        slope_numer = -two *sl2_numer;
        slope_denom = sl2_denom;
      }
      else {
        slope_numer = two *sl2_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // The third-order derivative is given by:
    //
    //              (2t*alpha - t*beta) * gamma
    //   y''' = -6 ------------------------------
    //                    beta^2 * delta
    //
    const Algebraic sl3_numer =
      (two * m_alg_r * sl_numer - m_alg_t * sl_denom) * sl2_numer;
    const Algebraic sl3_denom = sl_denom*sl_denom * sl2_denom;

    if (i == 3) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign(sl_denom) != NEGATIVE) {
        slope_numer = -6 * sl3_numer;
        slope_denom = sl3_denom;
      }
      else {
        slope_numer = 6 * sl3_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // \todo Handle higher-order derivatives as well.
    CGAL_error();
  }

private:
  /*! Compute the multiplicity of an intersection point.
   * \param arc The arc to intersect with.
   * \param p The intersection point.
   * \return The multiplicity of the intersection point.
   */
  unsigned int multiplicity_of_intersection_point(const Self& xcv,
                                                  const Point_2& p) const {
    CGAL_assertion(! is_special_segment() || ! xcv.is_special_segment());

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic slope1_numer, slope1_denom;
    Algebraic slope2_numer, slope2_denom;

    derive_by_x_at(p, 1, slope1_numer, slope1_denom);
    xcv.derive_by_x_at(p, 1, slope2_numer, slope2_denom);

    if (CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom) !=
        EQUAL) {
      // Different slopes at p - the mutiplicity of p is 1:
      return 1;
    }

    if (CGAL::sign(slope1_denom) != ZERO &&
        CGAL::sign(slope2_denom) != ZERO) {
      // The curves do not have a vertical slope at p.
      // Compare their second-order derivative by x:
      derive_by_x_at(p, 2, slope1_numer, slope1_denom);
      xcv.derive_by_x_at(p, 2, slope2_numer, slope2_denom);
    }
    else {
      // Both curves have a vertical slope at p.
      // Compare their second-order derivative by y:
      derive_by_y_at(p, 2, slope1_numer, slope1_denom);
      xcv.derive_by_y_at(p, 2, slope2_numer, slope2_denom);
    }

    if (CGAL::compare(slope1_numer*slope2_denom,
                      slope2_numer*slope1_denom) != EQUAL)
    {
      // Different curvatures at p - the mutiplicity of p is 2:
      return 2;
    }

    // If we reached here, the multiplicity of the intersection point is 3:
    return 3;
  }
  //@}

};

/*! Exporter for x-monotone conic arcs.
 */
template <typename Conic_arc_2>
std::ostream& operator<<(std::ostream& os,
                         const _Conic_x_monotone_arc_2<Conic_arc_2>& xcv)
{
  // Output the supporting conic curve.
  os << "{" << CGAL::to_double(xcv.r()) << "*x^2 + "
     << CGAL::to_double(xcv.s()) << "*y^2 + "
     << CGAL::to_double(xcv.t()) << "*xy + "
     << CGAL::to_double(xcv.u()) << "*x + "
     << CGAL::to_double(xcv.v()) << "*y + "
     << CGAL::to_double(xcv.w()) << "}";

  // Output the endpoints.
  os << " : (" << CGAL::to_double(xcv.source().x()) << ","
     << CGAL::to_double(xcv.source().y()) << ") ";

  if (xcv.orientation() == CLOCKWISE) os << "--cw-->";
  else if (xcv.orientation() == COUNTERCLOCKWISE) os << "--ccw-->";
  else os << "--l-->";

  os << " (" << CGAL::to_double(xcv.target().x()) << ","
     << CGAL::to_double(xcv.target().y()) << ")";

  return os;
}

} //namespace CGAL

#endif
