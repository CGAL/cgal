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
class _Conic_x_monotone_arc_2 : private ConicArc {
public:
  typedef ConicArc                                Conic_arc_2;
  typedef _Conic_x_monotone_arc_2<Conic_arc_2>    Self;

  typedef typename Conic_arc_2::Alg_kernel        Alg_kernel;
  typedef typename Conic_arc_2::Algebraic         Algebraic;

  typedef typename Conic_arc_2::Point_2           Point_2;
  typedef typename Conic_arc_2::Conic_point_2     Conic_point_2;

  // Type definition for the intersection points mapping.
  typedef typename Conic_point_2::Conic_id        Conic_id;
  typedef std::pair<Conic_id, Conic_id>           Conic_pair;
  typedef std::pair<Conic_point_2, unsigned int>  Intersection_point;
  typedef std::list<Intersection_point>           Intersection_list;

  using Conic_arc_2::sign_of_extra_data;
  using Conic_arc_2::is_between_endpoints;
  using Conic_arc_2::is_strictly_between_endpoints;
  using Conic_arc_2::conic_get_y_coordinates;

  using Conic_arc_2::IS_VALID;
  using Conic_arc_2::IS_FULL_CONIC;

  using Conic_arc_2::flag_mask;
  using Conic_arc_2::reset_flags;
  using Conic_arc_2::set_flag;
  using Conic_arc_2::reset_flag;
  using Conic_arc_2::flip_flag;
  using Conic_arc_2::test_flag;

  /*! \struct Less functor for Conic_pair.
   */
  struct Less_conic_pair {
    bool operator()(const Conic_pair& cp1, const Conic_pair& cp2) const {
      // Compare the pairs of IDs lexicographically.
      return ((cp1.first < cp2.first) ||
              ((cp1.first == cp2.first) && (cp1.second < cp2.second)));
    }
  };

  typedef std::map<Conic_pair, Intersection_list, Less_conic_pair>
                                                  Intersection_map;
  typedef typename Intersection_map::value_type   Intersection_map_entry;
  typedef typename Intersection_map::iterator     Intersection_map_iterator;

protected:
  typedef Conic_arc_2                             Base;

  typedef typename Conic_arc_2::Integer           Integer;
  typedef typename Conic_arc_2::Nt_traits         Nt_traits;
  typedef typename Conic_arc_2::Rat_kernel        Rat_kernel;

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

  Algebraic alg_r;      // The coefficients of the supporting conic curve:
  Algebraic alg_s;      //
  Algebraic alg_t;      //   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0 ,
  Algebraic alg_u;      //
  Algebraic alg_v;      // converted to algebraic numbers.
  Algebraic alg_w;      //

  Conic_id m_id;         // The ID number of the supporting conic curve.

public:
  /// \name Constrcution methods.
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
    alg_r(arc.alg_r),
    alg_s(arc.alg_s),
    alg_t(arc.alg_t),
    alg_u(arc.alg_u),
    alg_v(arc.alg_v),
    alg_w(arc.alg_w),
    m_id(arc.m_id)
  {}

  /*! Construct an x-monotone arc from a conic arc.
   * \param arc The given (base) arc.
   * \pre The given arc is x-monotone.
   */
  _Conic_x_monotone_arc_2(const Base& arc) :
    Base(arc),
    m_id()
  {
    CGAL_precondition(arc.is_valid() && arc.is_x_monotone());
    set();
  }

  /*! Construct an x-monotone arc from a conic arc.
   * \param arc The given (base) arc.
   * \param id The ID of the base arc.
   */
  _Conic_x_monotone_arc_2(const Base& arc, const Conic_id& id) :
    Base(arc),
    m_id(id)
  {
    CGAL_precondition(arc.is_valid() && id.is_valid());
    set();
  }

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
  {
    CGAL_precondition(arc.is_valid() && id.is_valid());

    // Set the two endpoints.
    this->m_source = source;
    this->m_target = target;

    set();
  }

  /*! Construct a special segment connecting to given endpoints (for the usage
   * of the landmarks point-location strategy).
   * \param source The source point.
   * \param target The target point.
   */
  _Conic_x_monotone_arc_2(const Point_2& source, const Point_2& target) :
    Base(source, target)
  {
    set_flag(DEGREE_1);

    Alg_kernel ker;
    auto cmp_xy = ker.compare_xy_2_object();
    Comparison_result dir_res = cmp_xy(this->m_source, this->m_target);
    if (dir_res == SMALLER) set_flag(IS_DIRECTED_RIGHT);

    // Check if the segment is vertical.
    if (CGAL::sign(this->m_extra_data->b) == ZERO)
      set_flag(IS_VERTICAL_SEGMENT);

    // Mark that this is a special segment.
    set_flag(IS_SPECIAL_SEGMENT);
  }

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
  {
    // Make sure the two endpoints lie on the supporting line.
    CGAL_precondition
      (CGAL::sign(a * source.x() + b * source.y() + c) == CGAL::ZERO);

    CGAL_precondition
      (CGAL::sign(a * target.x() + b * target.y() + c) == CGAL::ZERO);

    // Set the basic properties and clear the m_info bits.
    this->m_source = source;
    this->m_target = target;
    this->m_orient = COLLINEAR;
    reset_flags();              // inavlid arc

    // Check if the arc is directed right (the target is lexicographically
    // greater than the source point), or to the left.
    Alg_kernel ker;
    Comparison_result res =
      ker.compare_x_2_object()(this->m_source, this->m_target);

    set_flag(IS_VALID);
    set_flag(DEGREE_1);
    if (res == EQUAL) {
      // Mark that the segment is vertical.
      set_flag(IS_VERTICAL_SEGMENT);

      // Compare the endpoints lexicographically.
      res = ker.compare_y_2_object()(this->m_source, this->m_target);

      CGAL_precondition (res != EQUAL);
      if (res == EQUAL) {
        reset_flags();          // inavlid arc
        return;
      }
    }

    if (res == SMALLER) set_flag(IS_DIRECTED_RIGHT);

    // Store the coefficients of the line.
    this->m_extra_data = new typename Base::Extra_data;
    this->m_extra_data->a = a;
    this->m_extra_data->b = b;
    this->m_extra_data->c = c;
    this->m_extra_data->side = ZERO;

    // Mark that this is a special segment.
    set_flag(IS_SPECIAL_SEGMENT);
  }

  /*! Assignment operator.
   * \param arc The copied arc.
   */
  const Self& operator=(const Self& arc) {
    CGAL_precondition (arc.is_valid());

    if (this == &arc) return (*this);

    // Copy the base arc.
    Base::operator= (arc);

    // Set the rest of the properties.
    alg_r = arc.alg_r;
    alg_s = arc.alg_s;
    alg_t = arc.alg_t;
    alg_u = arc.alg_u;
    alg_v = arc.alg_v;
    alg_w = arc.alg_w;

    m_id = arc.m_id;

    return (*this);
  }
  //@}

  /// \name Accessing the arc properties.
  //@{

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

  /*! Return true iff the conic arc is directed right iexicographically.
   */
  bool is_directed_right() const
  { return test_flag(IS_DIRECTED_RIGHT); }

  /*! Obtain a bounding box for the conic arc.
   * \return The bounding box.
   */
  Bbox_2 bbox() const { return Base::bbox(); }
  //@}

  /// \name Predicates.
  //@{

  /*! Check if the conic arc is a vertical segment.
   */
  bool is_vertical() const
  { return test_flag(IS_VERTICAL_SEGMENT); }

  /*! Check whether the given point lies on the arc.
   * \param p The qury point.
   * \param (true) if p lies on the arc; (false) otherwise.
   */
  bool contains_point(const Conic_point_2& p) const {
    // First check if p lies on the supporting conic. We first check whether
    // it is one of p's generating conic curves.
    bool p_on_conic(false);

    if (p.is_generating_conic(m_id)) p_on_conic = true;
    else {
      // Check whether p satisfies the supporting conic equation.
      p_on_conic = is_on_supporting_conic(p.x(), p.y());

      if (p_on_conic) {
        // As p lies on the supporting conic of our arc, add its ID to
        // the list of generating conics for p.
        Conic_point_2& p_non_const = const_cast<Conic_point_2&>(p);
        p_non_const.set_generating_conic(m_id);
      }
    }

    if (! p_on_conic) return false;

    // Check if p is between the endpoints of the arc.
    return is_between_endpoints(p);
  }
  //@}

  /// \name Constructing points on the arc.
  //@{

  /*! Compute a point on the arc with the same x-coordiante as the given point.
   * \param p The given point.
   * \pre The arc is not vertical and p is in the x-range of the arc.
   * \return A point on the arc with the same x-coordiante as p.
   */
  Point_2 point_at_x(const Point_2& p) const {
    // Make sure that p is in the x-range of the arc.
    CGAL_precondition(! test_flag(IS_VERTICAL_SEGMENT));

    CGAL_precondition_code(Alg_kernel ker;);

    CGAL_precondition(ker.compare_x_2_object()(p, left()) != SMALLER &&
                      ker.compare_x_2_object()(p, right()) != LARGER);

    if (is_special_segment()) {
      // In case of a special segment, the equation of the supported line
      // (a*x + b*y + c) = 0 is stored with the extra data field, and we
      // simply have:
      Algebraic _y = -(this->m_extra_data->a*p.x() + this->m_extra_data->c) /
        this->m_extra_data->b;

      // Return the computed point.
      return Point_2(p.x(), _y);
    }

    // Compute the y-coordinate according to the degree of the supporting
    // conic curve.
    Nt_traits nt_traits;
    Algebraic y;

    if ((this->m_info & DEGREE_MASK) == flag_mask(DEGREE_1)) {
      // In case of a linear curve, the y-coordinate is a simple linear
      // expression of x(p) (note that v is not 0 as the arc is not vertical):
      //   y = -(u*x(p) + w) / v
      y = -(alg_u*p.x() + alg_w) / alg_v;
    }
    else if (this->m_orient == COLLINEAR) {
      CGAL_assertion (this->m_extra_data != nullptr);

      // In this case the equation of the supporting line is given by the
      // extra data structure.
      y = -(this->m_extra_data->a * p.x() +
            this->m_extra_data->c) / this->m_extra_data->b;
    }
    else {
      CGAL_assertion((this->m_info & DEGREE_MASK) == flag_mask(DEGREE_2));

      // In this case the y-coordinate is one of solutions to the quadratic
      // equation:
      //  s*y^2 + (t*x(p) + v)*y + (r*x(p)^2 + u*x(p) + w) = 0
      Algebraic A = alg_s;
      Algebraic B = alg_t*p.x() + alg_v;
      Algebraic C = (alg_r*p.x() + alg_u)*p.x() + alg_w;

      if (CGAL::sign(this->m_s) == ZERO) {
        // In this case A is 0 and we have a linear equation.
        CGAL_assertion (CGAL::sign (B) != ZERO);

        y = -C / B;
      }
      else {
        // Solve the quadratic equation.
        Algebraic  disc = B*B - 4*A*C;

        CGAL_assertion (CGAL::sign (disc) != NEGATIVE);

        // We take either the root involving -sqrt(disc) or +sqrt(disc)
        // based on the information flags.
        y = (test_flag(PLUS_SQRT_DISC_ROOT)) ?
          (nt_traits.sqrt (disc) - B) / (2*A) :
          -(B + nt_traits.sqrt (disc)) / (2*A);
      }
    }

    // Return the computed point.
    return Point_2(p.x(), y);
  }

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

  /*! Compare to arcs immediately to the right of their intersection point.
   * \param arc The compared arc.
   * \param p The reference intersection point.
   * \return The relative position of the arcs to the right of p.
   * \pre Both arcs we compare are not vertical segments.
   */
  Comparison_result compare_to_right(const Self& arc, const Conic_point_2& p)
    const
  {
    CGAL_precondition(! test_flag(IS_VERTICAL_SEGMENT) &&
                      ! arc.test_flag(IS_VERTICAL_SEGMENT));

    // In case one arc is facing upwards and another facing downwards, it is
    // clear that the one facing upward is above the one facing downwards.
    if (has_same_supporting_conic(arc)) {
      if (test_flag(FACING_UP) && arc.test_flag(FACING_DOWN))
        return LARGER;
      else if (test_flag(FACING_DOWN) && arc.test_flag(FACING_UP))
        return SMALLER;

      // In this case the two arcs overlap.
      CGAL_assertion((this->m_info & FACING_MASK) == (arc.m_info & FACING_MASK));
      return EQUAL;
    }

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic slope1_numer, slope1_denom;
    Algebraic slope2_numer, slope2_denom;

    derive_by_x_at(p, 1, slope1_numer, slope1_denom);
    arc.derive_by_x_at(p, 1, slope2_numer, slope2_denom);

    // Check if any of the slopes is vertical.
    const bool is_vertical_slope1 = (CGAL::sign(slope1_denom) == ZERO);
    const bool is_vertical_slope2 = (CGAL::sign(slope2_denom) == ZERO);

    if (! is_vertical_slope1 && !is_vertical_slope2) {
      // The two derivatives at p are well-defined: use them to determine
      // which arc is above the other (the one with a larger slope is below).
      Comparison_result slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the second-order derivative.
      derive_by_x_at(p, 2, slope1_numer, slope1_denom);
      arc.derive_by_x_at(p, 2, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the third-order derivative.
      derive_by_x_at(p, 3, slope1_numer, slope1_denom);
      arc.derive_by_x_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      // \todo Handle higher-order derivatives:
      CGAL_assertion(slope_res != EQUAL);

      return slope_res;
    }
    else if (! is_vertical_slope2) {
      // The first arc has a vertical slope at p: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion((this->m_info & FACING_MASK) != 0);

      return (test_flag(FACING_UP)) ? LARGER : SMALLER;
    }
    else if (! is_vertical_slope1) {
      // The second arc has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion((arc.m_info & FACING_MASK) != 0);

      return (arc.test_flag(FACING_UP)) ? SMALLER : LARGER;
    }

    // The two arcs have vertical slopes at p_int:
    // First check whether one is facing up and one down. In this case the
    // comparison result is trivial.
    if (test_flag(FACING_UP) && arc.test_flag(FACING_DOWN))
      return LARGER;
    else if (test_flag(FACING_DOWN) && arc.test_flag(FACING_UP))
      return SMALLER;

    // Compute the second-order derivative by y and act according to it.
    derive_by_y_at(p, 2, slope1_numer, slope1_denom);
    arc.derive_by_y_at(p, 2, slope2_numer, slope2_denom);

    Comparison_result slope_res =
      CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

    // If necessary, use the third-order derivative by y.
    if (slope_res == EQUAL) {
      // \todo Check this!
      derive_by_y_at(p, 3, slope1_numer, slope1_denom);
      arc.derive_by_y_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);
    }

    // \todo Handle higher-order derivatives:
    CGAL_assertion(slope_res != EQUAL);

    if (test_flag(FACING_UP) && arc.test_flag(FACING_UP)) {
      // Both are facing up.
      return ((slope_res == LARGER) ? SMALLER : LARGER);
    }
    // Both are facing down.
    return slope_res;
  }

  /*! Compare to arcs immediately to the leftt of their intersection point.
   * \param arc The compared arc.
   * \param p The reference intersection point.
   * \return The relative position of the arcs to the left of p.
   * \pre Both arcs we compare are not vertical segments.
   */
  Comparison_result compare_to_left(const Self& arc, const Conic_point_2& p)
    const
  {
    CGAL_precondition(! test_flag(IS_VERTICAL_SEGMENT) &&
                      ! arc.test_flag(IS_VERTICAL_SEGMENT));

    // In case one arc is facing upwards and another facing downwards, it is
    // clear that the one facing upward is above the one facing downwards.
    if (has_same_supporting_conic (arc)) {
      if (test_flag(FACING_UP) && arc.test_flag(FACING_DOWN))
        return LARGER;
      else if (test_flag(FACING_DOWN) && arc.test_flag(FACING_UP))
        return SMALLER;

      // In this case the two arcs overlap.
      CGAL_assertion((this->m_info & FACING_MASK) == (arc.m_info & FACING_MASK));

      return EQUAL;
    }

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic slope1_numer, slope1_denom;
    Algebraic slope2_numer, slope2_denom;

    derive_by_x_at(p, 1, slope1_numer, slope1_denom);
    arc.derive_by_x_at(p, 1, slope2_numer, slope2_denom);

    // Check if any of the slopes is vertical.
    const bool is_vertical_slope1 = (CGAL::sign (slope1_denom) == ZERO);
    const bool is_vertical_slope2 = (CGAL::sign (slope2_denom) == ZERO);

    if (! is_vertical_slope1 && !is_vertical_slope2) {
      // The two derivatives at p are well-defined: use them to determine
      // which arc is above the other (the one with a larger slope is below).
      Comparison_result slope_res = CGAL::compare(slope2_numer*slope1_denom,
                                                  slope1_numer*slope2_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the second-order derivative.
      derive_by_x_at(p, 2, slope1_numer, slope1_denom);
      arc.derive_by_x_at(p, 2, slope2_numer, slope2_denom);

      slope_res = CGAL::compare(slope1_numer*slope2_denom,
                                slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the third-order derivative.
      derive_by_x_at(p, 3, slope1_numer, slope1_denom);
      arc.derive_by_x_at(p, 3, slope2_numer, slope2_denom);

      slope_res = CGAL::compare(slope2_numer*slope1_denom,
                                slope1_numer*slope2_denom);

      // \todo Handle higher-order derivatives:
      CGAL_assertion(slope_res != EQUAL);

      return slope_res;
    }
    else if (! is_vertical_slope2) {
      // The first arc has a vertical slope at p: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion((this->m_info & FACING_MASK) != 0);

      return (test_flag(FACING_UP)) ? LARGER : SMALLER;
    }
    else if (! is_vertical_slope1) {
      // The second arc has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion((arc.m_info & FACING_MASK) != 0);

      return (arc.test_flag(FACING_UP)) ? SMALLER : LARGER;
    }

    // The two arcs have vertical slopes at p_int:
    // First check whether one is facing up and one down. In this case the
    // comparison result is trivial.
    if (test_flag(FACING_UP) && arc.test_flag(FACING_DOWN))
      return LARGER;
    else if (test_flag(FACING_DOWN) && arc.test_flag(FACING_UP))
      return SMALLER;

    // Compute the second-order derivative by y and act according to it.
    derive_by_y_at(p, 2, slope1_numer, slope1_denom);
    arc.derive_by_y_at(p, 2, slope2_numer, slope2_denom);

    Comparison_result slope_res =
      CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);

    // If necessary, use the third-order derivative by y.
    if (slope_res == EQUAL) {
      // \todo Check this!
      derive_by_y_at(p, 3, slope1_numer, slope1_denom);
      arc.derive_by_y_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);
    }

    // \todo Handle higher-order derivatives:
    CGAL_assertion(slope_res != EQUAL);

    // Check whether both are facing up.
    if (test_flag(FACING_UP) && arc.test_flag(FACING_UP))
      return ((slope_res == LARGER) ? SMALLER : LARGER);

    // Both are facing down.
    return slope_res;
  }

  /*! Compute the intersections with the given arc.
   * \param arc The given intersecting arc.
   * \param inter_map Maps conic pairs to lists of their intersection points.
   * \param oi The output iterator.
   * \return The past-the-end iterator.
   */
  template <typename OutputIterator>
  OutputIterator intersect(const Self& arc, Intersection_map& inter_map,
                           OutputIterator oi) const
  {
    typedef boost::variant<Intersection_point, Self>  Intersection_result;

    if (has_same_supporting_conic(arc)) {
      // Check for overlaps between the two arcs.
      Self overlap;

      if (compute_overlap(arc, overlap)) {
        // There can be just a single overlap between two x-monotone arcs:
        *oi++ = Intersection_result(overlap);
        return oi;
      }

      // In case there is not overlap and the supporting conics are the same,
      // there cannot be any intersection points, unless the two arcs share
      // an end point.
      // Note that in this case we do not define the multiplicity of the
      // intersection points we report.
      Alg_kernel ker;

      if (ker.equal_2_object()(left(), arc.left())) {
        Intersection_point ip(left(), 0);
        *oi++ = Intersection_result(ip);
      }

      if (ker.equal_2_object()(right(), arc.right())) {
        Intersection_point ip(right(), 0);
        *oi++ = Intersection_result(ip);
      }

      return oi;
    }

    // Search for the pair of supporting conics in the map (the first conic
    // ID in the pair should be smaller than the second one, to guarantee
    // uniqueness).
    Conic_pair conic_pair;
    Intersection_map_iterator map_iter;
    Intersection_list inter_list;
    bool invalid_ids = false;

    if (m_id.is_valid() && arc.m_id.is_valid()) {
      if (m_id < arc.m_id) conic_pair = Conic_pair(m_id, arc.m_id);
      else conic_pair = Conic_pair(arc.m_id, m_id);
      map_iter = inter_map.find(conic_pair);
    }
    else {
      // In case one of the IDs is invalid, we do not look in the map neither
      // we cache the results.
      map_iter = inter_map.end();
      invalid_ids = true;
    }

    if (map_iter == inter_map.end()) {
      // In case the intersection points between the supporting conics have
      // not been computed before, compute them now and store them in the map.
      intersect_supporting_conics(arc, inter_list);

      if (! invalid_ids) inter_map[conic_pair] = inter_list;
    }
    else {
      // Obtain the precomputed intersection points from the map.
      inter_list = (*map_iter).second;
    }

    // Go over the list of intersection points and report those that lie on
    // both x-monotone arcs.
    for (auto iter = inter_list.begin(); iter != inter_list.end(); ++iter) {
      if (is_between_endpoints((*iter).first) &&
          arc.is_between_endpoints((*iter).first))
      {
        *oi++ = Intersection_result(*iter);
      }
    }

    return oi;
  }
  //@}

  /// \name Constructing x-monotone arcs.
  //@{

  /*! Split the arc into two at a given split point.
   * \param p The split point.
   * \param c1 Output: The first resulting arc, lying to the left of p.
   * \param c2 Output: The first resulting arc, lying to the right of p.
   * \pre p lies in the interior of the arc (not one of its endpoints).
   */
  void split(const Conic_point_2& p, Self& c1, Self& c2) const {
    // Make sure that p lies on the interior of the arc.
    CGAL_precondition_code(Alg_kernel ker);

    CGAL_precondition(this->contains_point (p) &&
                      ! ker.equal_2_object()(p, this->m_source) &&
                      ! ker.equal_2_object()(p, this->m_target));

    // Make copies of the current arc.
    c1 = *this;
    c2 = *this;

    // Assign the endpoints of the arc.
    if (test_flag(IS_DIRECTED_RIGHT)) {
      // The arc is directed from left to right, so p becomes c1's target
      // and c2's source.
      c1.m_target = p;
      c2.m_source = p;

      if (! p.is_generating_conic(m_id)) {
        c1.m_target.set_generating_conic(m_id);
        c2.m_source.set_generating_conic(m_id);
      }
    }
    else {
      // The arc is directed from right to left, so p becomes c2's target
      // and c1's source.
      c1.m_source = p;
      c2.m_target = p;

      if (! p.is_generating_conic(m_id)) {
        c1.m_source.set_generating_conic(m_id);
        c2.m_target.set_generating_conic(m_id);
      }
    }
  }

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

  /*! Trim the arc given its new endpoints.
   * \param ps The new source point.
   * \param pt The new target point.
   * \return The new trimmed arc.
   * \pre Both ps and pt lies on the arc and must conform with the current
   *      direction of the arc.
   */
  Self trim(const Conic_point_2& ps, const Conic_point_2& pt) const {
    // Make sure that both ps and pt lie on the arc.
    CGAL_precondition(this->contains_point(ps) && this->contains_point(pt));

    // Make sure that the endpoints conform with the direction of the arc.
    Self arc = *this;
    Alg_kernel ker;

    if (! ((test_flag(IS_DIRECTED_RIGHT) &&
            (ker.compare_xy_2_object()(ps, pt) == SMALLER)) ||
           (test_flag(IS_DIRECTED_RIGHT) &&
            (ker.compare_xy_2_object()(ps, pt) == LARGER))))
    {
      // We are allowed to change the direction only in case of a segment.
      CGAL_assertion(this->m_orient == COLLINEAR);
      arc.flip_flag(IS_DIRECTED_RIGHT);
    }

    // Make a copy of the current arc and assign its endpoints.
    if (! ker.equal_2_object()(ps, this->m_source)) {
      arc.m_source = ps;

      if (! ps.is_generating_conic(m_id))
        arc.m_source.set_generating_conic(m_id);
    }

    if (! ker.equal_2_object()(pt, this->m_target)) {
      arc.m_target = pt;

      if (! pt.is_generating_conic(m_id))
        arc.m_target.set_generating_conic(m_id);
    }

    return arc;
  }

  /*! Check whether the two arcs are equal (have the same graph).
   * \param arc The compared arc.
   * \return (true) if the two arcs have the same graph; (false) otherwise.
   */
  bool equals(const Self& arc) const {
    // The two arc must have the same supporting conic curves.
    if (! has_same_supporting_conic (arc)) return false;

    // Check that the arc endpoints are the same.
    Alg_kernel ker;

    if (this->m_orient == COLLINEAR) {
      CGAL_assertion(arc.m_orient == COLLINEAR);
      return((ker.equal_2_object()(this->m_source, arc.m_source) &&
              ker.equal_2_object()(this->m_target, arc.m_target)) ||
              (ker.equal_2_object()(this->m_source, arc.m_target) &&
               ker.equal_2_object()(this->m_target, arc.m_source)));
    }

    if (this->m_orient == arc.m_orient) {
      // Same orientation - the source and target points must be the same.
      return (ker.equal_2_object()(this->m_source, arc.m_source) &&
              ker.equal_2_object()(this->m_target, arc.m_target));
    }
    else {
      // Reverse orientation - the source and target points must be swapped.
      return (ker.equal_2_object()(this->m_source, arc.m_target) &&
              ker.equal_2_object()(this->m_target, arc.m_source));
    }
  }

  /*! Check whether it is possible to merge the arc with the given arc.
   * \param arc The query arc.
   * \return (true) if it is possible to merge the two arcs;
   *         (false) otherwise.
   */
  bool can_merge_with(const Self& arc) const {
    // In order to merge the two arcs, they should have the same supporting
    // conic.
    if (! has_same_supporting_conic(arc)) return false;

    // Check if the left endpoint of one curve is the right endpoint of the
    // other.
    Alg_kernel ker;
    return (ker.equal_2_object() (right(), arc.left()) ||
            ker.equal_2_object() (left(), arc.right()));
  }

  /*! Merge the current arc with the given arc.
   * \param arc The arc to merge with.
   * \pre The two arcs are mergeable.
   */
  void merge(const Self& arc) {
    CGAL_precondition (this->can_merge_with (arc));

    // Check if we should extend the arc to the left or to the right.
    Alg_kernel ker;

    if (ker.equal_2_object() (right(), arc.left())) {
      // Extend the arc to the right.
      if (test_flag(IS_DIRECTED_RIGHT)) this->m_target = arc.right();
      else this->m_source = arc.right();
    }
    else {
      CGAL_precondition(ker.equal_2_object() (left(), arc.right()));

      // Extend the arc to the left.
      if (test_flag(IS_DIRECTED_RIGHT))
        this->m_source = arc.left();
      else
        this->m_target = arc.left();
    }
  }

  bool is_upper() const { return test_flag(FACING_UP); }

  bool is_lower() const { return test_flag(FACING_DOWN); }
  //@}

private:
  /// \name Auxiliary (private) functions.
  //@{

  /*! Set the properties of the x-monotone conic arc (for the usage of the
   * constructors).
   */
  void set() {
    // Convert the coefficients of the supporting conic to algebraic numbers.
    Nt_traits nt_traits;

    alg_r = nt_traits.convert(this->m_r);
    alg_s = nt_traits.convert(this->m_s);
    alg_t = nt_traits.convert(this->m_t);
    alg_u = nt_traits.convert(this->m_u);
    alg_v = nt_traits.convert(this->m_v);
    alg_w = nt_traits.convert(this->m_w);

    // Set the generating conic ID for the source and target points.
    this->m_source.set_generating_conic(m_id);
    this->m_target.set_generating_conic(m_id);

    // Update the m_info bits.
    set_flag(IS_VALID);
    reset_flag(IS_FULL_CONIC);

    // Check if the arc is directed right (the target is lexicographically
    // greater than the source point), or to the left.
    Alg_kernel ker;
    Comparison_result dir_res =
      ker.compare_xy_2_object()(this->m_source, this->m_target);

    CGAL_assertion(dir_res != EQUAL);

    if (dir_res == SMALLER) set_flag(IS_DIRECTED_RIGHT);

    // Compute the degree of the underlying conic.
    if ((CGAL::sign(this->m_r) != ZERO) ||
        (CGAL::sign(this->m_s) != ZERO) ||
        (CGAL::sign(this->m_t) != ZERO))
    {
      set_flag(DEGREE_2);

      if (this->m_orient == COLLINEAR) {
        set_flag(IS_SPECIAL_SEGMENT);

        // Check whether the arc is a vertical segment:
        if (ker.compare_x_2_object()(this->m_source, this->m_target) == EQUAL)
          set_flag(IS_VERTICAL_SEGMENT);

        return;
      }
    }
    else {
      CGAL_assertion(CGAL::sign(this->m_u) != ZERO ||
                     CGAL::sign(this->m_v) != ZERO);

      if (CGAL::sign(this->m_v) == ZERO) {

        // The supporting curve is of the form: _u*x + _w = 0
        set_flag(IS_VERTICAL_SEGMENT);
      }

      set_flag(DEGREE_1);

      return;
    }

    if (this->m_orient == COLLINEAR) return;

    // Compute a midpoint between the source and the target and get the y-value
    // of the arc at its x-coordiante.
    Point_2 p_mid =
      ker.construct_midpoint_2_object()(this->m_source, this->m_target);
    Algebraic ys[2];
    CGAL_assertion_code(int n_ys = )
      conic_get_y_coordinates(p_mid.x(), ys);

    CGAL_assertion(n_ys != 0);

    // Check which solution lies on the x-monotone arc.
    Point_2 p_arc_mid(p_mid.x(), ys[0]);

    if (is_strictly_between_endpoints(p_arc_mid)) {
      // Mark that we should use the -sqrt(disc) root for points on this
      // x-monotone arc.
      reset_flag(PLUS_SQRT_DISC_ROOT);
    }
    else {
      CGAL_assertion (n_ys == 2);
      p_arc_mid = Point_2 (p_mid.x(), ys[1]);

      CGAL_assertion(is_strictly_between_endpoints(p_arc_mid));

      // Mark that we should use the +sqrt(disc) root for points on this
      // x-monotone arc.
      set_flag(PLUS_SQRT_DISC_ROOT);
    }

    // Check whether the conic is facing up or facing down:
    // Check whether the arc (which is x-monotone of degree 2) lies above or
    // below the segement that contects its two end-points (x1,y1) and (x2,y2).
    // To do that, we find the y coordinate of a point on the arc whose x
    // coordinate is (x1+x2)/2 and compare it to (y1+y2)/2.
    Comparison_result res = ker.compare_y_2_object() (p_arc_mid, p_mid);

    // If the arc is above the connecting segment, so it is facing upwards.
    if (res == LARGER) set_flag(FACING_UP);
    // If the arc is below the connecting segment, so it is facing downwards.
    else if (res == SMALLER) set_flag(FACING_DOWN);
  }

  /*! Check if the arc is a special segment connecting two algebraic endpoints
   * (and has no undelying integer conic coefficients).
   */
  bool is_special_segment() const
  { return test_flag(IS_SPECIAL_SEGMENT); }

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
      _sign = CGAL::sign((alg_r*px + alg_t*py + alg_u) * px +
                         (alg_s*py + alg_v) * py + alg_w);
    }
    else {
      // Check whether p satisfies the equation of the line stored with the
      // extra data.
      _sign = sign_of_extra_data(px, py);
    }

    return (_sign == ZERO);
  }

  /*! Check whether the two arcs have the same supporting conic.
   * \param arc The compared arc.
   * \return (true) if the two supporting conics are the same.
   */
  bool has_same_supporting_conic(const Self& arc) const {
    // Check if the two arcs originate from the same conic:
    if (m_id == arc.m_id && m_id.is_valid() && arc.m_id.is_valid()) return true;

    // In case both arcs are collinear, check if they have the same
    // supporting lines.
    if ((this->m_orient == COLLINEAR) && (arc.m_orient == COLLINEAR)) {
      // Construct the two supporting lines and compare them.
      Alg_kernel ker;
      auto construct_line = ker.construct_line_2_object();
      typename Alg_kernel::Line_2 l1 =
        construct_line(this->m_source, this->m_target);
      typename Alg_kernel::Line_2 l2 =
        construct_line(arc.m_source, arc.m_target);
      auto equal = ker.equal_2_object();

      if (equal(l1, l2)) return true;

      // Try to compare l1 with the opposite of l2.
      l2 = construct_line(arc.m_target, arc.m_source);

      return equal(l1, l2);
    }
    else if ((this->m_orient == COLLINEAR) || (arc.m_orient == COLLINEAR)) {
      // Only one arc is collinear, so the supporting curves cannot be the
      // same:
      return false;
    }

    // Check whether the coefficients of the two supporting conics are equal
    // up to a constant factor.
    Integer factor1 = 1;
    Integer factor2 = 1;

    if (CGAL::sign(this->m_r) != ZERO) factor1 = this->m_r;
    else if (CGAL::sign(this->m_s) != ZERO) factor1 = this->m_s;
    else if (CGAL::sign(this->m_t) != ZERO) factor1 = this->m_t;
    else if (CGAL::sign(this->m_u) != ZERO) factor1 = this->m_u;
    else if (CGAL::sign(this->m_v) != ZERO) factor1 = this->m_v;
    else if (CGAL::sign(this->m_w) != ZERO) factor1 = this->m_w;

    if (CGAL::sign(arc.m_r) != ZERO) factor2 = arc.m_r;
    else if (CGAL::sign(arc.m_s) != ZERO) factor2 = arc.m_s;
    else if (CGAL::sign(arc.m_t) != ZERO) factor2 = arc.m_t;
    else if (CGAL::sign(arc.m_u) != ZERO) factor2 = arc.m_u;
    else if (CGAL::sign(arc.m_v) != ZERO) factor2 = arc.m_v;
    else if (CGAL::sign(arc.m_w) != ZERO) factor2 = arc.m_w;

    return (CGAL::compare(this->m_r * factor2, arc.m_r * factor1) == EQUAL &&
            CGAL::compare(this->m_s * factor2, arc.m_s * factor1) == EQUAL &&
            CGAL::compare(this->m_t * factor2, arc.m_t * factor1) == EQUAL &&
            CGAL::compare(this->m_u * factor2, arc.m_u * factor1) == EQUAL &&
            CGAL::compare(this->m_v * factor2, arc.m_v * factor1) == EQUAL &&
            CGAL::compare(this->m_w * factor2, arc.m_w * factor1) == EQUAL);
  }

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
    const Algebraic sl_numer = two*alg_r*p.x() + alg_t*p.y() + alg_u;
    const Algebraic sl_denom = two*alg_s*p.y() + alg_t*p.x() + alg_v;

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
    const Algebraic sl2_numer = alg_s * sl_numer*sl_numer -
      alg_t * sl_numer*sl_denom + alg_r * sl_denom*sl_denom;
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
      (two * alg_s * sl_numer - alg_t * sl_denom) * sl2_numer;
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
        if (CGAL::sign (this->m_extra_data->a) != NEGATIVE) {
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
    const Algebraic sl_numer = two*alg_s*p.y() + alg_t*p.x() + alg_v;
    const Algebraic sl_denom = two*alg_r*p.x() + alg_t*p.y() + alg_u;

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
    const Algebraic sl2_numer = alg_r * sl_numer*sl_numer -
      alg_t * sl_numer*sl_denom + alg_s * sl_denom*sl_denom;
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
      (two * alg_r * sl_numer - alg_t * sl_denom) * sl2_numer;
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

  /*! Compute the overlap with a given arc, which is supposed to have the same
   * supporting conic curve as this arc.
   * \param arc The given arc.
   * \param overlap Output: The overlapping arc (if any).
   * \return Whether we found an overlap.
   */
  bool compute_overlap(const Self& arc, Self& overlap) const {
    std::cout << "compute_overlap()\n";
    // Check if the two arcs are identical.
    if (equals(arc)) {
      overlap = arc;
      return true;
    }

    if (is_strictly_between_endpoints(arc.left())) {
      if (is_strictly_between_endpoints(arc.right())) {
        // Case 1 - *this:     +----------->
        //            arc:       +=====>
        overlap = arc;
        return true;
      }
      else {
        // Case 2 - *this:     +----------->
        //            arc:               +=====>
        overlap = *this;

        if (overlap.test_flag(IS_DIRECTED_RIGHT))
          overlap.m_source = arc.left();
        else
          overlap.m_target = arc.left();

        return true;
      }
    }
    else if (is_strictly_between_endpoints(arc.right())) {
      // Case 3 - *this:     +----------->
      //            arc:   +=====>
      overlap = *this;

      if (overlap.test_flag(IS_DIRECTED_RIGHT))
        overlap.m_target = arc.right();
      else
        overlap.m_source = arc.right();

      return true;
    }
    else if (arc.is_between_endpoints(this->m_source) &&
             arc.is_between_endpoints(this->m_target) &&
             (arc.is_strictly_between_endpoints(this->m_source) ||
              arc.is_strictly_between_endpoints(this->m_target)))
    {
      // Case 4 - *this:     +----------->
      //            arc:   +================>
      overlap = *this;
      return true;
    }

    // If we reached here, there are no overlaps:
    return false;
  }

  /*! Intersect the supporing conic curves of this arc and the given arc.
   * \param arc The arc to intersect with.
   * \param inter_list The list of intersection points.
   */
  void intersect_supporting_conics(const Self& arc,
                                   Intersection_list& inter_list) const {
    if (is_special_segment() && ! arc.is_special_segment()) {
      // If one of the arcs is a special segment, make sure it is (arc).
      arc.intersect_supporting_conics(*this, inter_list);
      return;
    }

    const int deg1 = ((this->m_info & DEGREE_MASK) == flag_mask(DEGREE_1)) ? 1 : 2;
    const int deg2 = ((arc.m_info & DEGREE_MASK) == flag_mask(DEGREE_1)) ? 1 : 2;
    Nt_traits nt_traits;
    Algebraic xs[4];
    int n_xs = 0;
    Algebraic ys[4];
    int n_ys = 0;

    if (arc.is_special_segment()) {
      // The second arc is a special segment (a*x + b*y + c = 0).
      if (is_special_segment()) {
        // Both arc are sepcial segment, so they have at most one intersection
        // point.
        Algebraic denom = this->m_extra_data->a * arc.m_extra_data->b -
          this->m_extra_data->b * arc.m_extra_data->a;

        if (CGAL::sign (denom) != CGAL::ZERO) {
          xs[0] = (this->m_extra_data->b * arc.m_extra_data->c -
                   this->m_extra_data->c * arc.m_extra_data->b) / denom;
          n_xs = 1;

          ys[0] = (this->m_extra_data->c * arc.m_extra_data->a -
                   this->m_extra_data->a * arc.m_extra_data->c) / denom;
          n_ys = 1;
        }
      }
      else {
        // Compute the x-coordinates of the intersection points.
        n_xs = compute_resultant_roots(nt_traits,
                                       alg_r, alg_s, alg_t,
                                       alg_u, alg_v, alg_w,
                                       deg1,
                                       arc.m_extra_data->a,
                                       arc.m_extra_data->b,
                                       arc.m_extra_data->c,
                                       xs);
        CGAL_assertion (n_xs <= 2);

        // Compute the y-coordinates of the intersection points.
        n_ys = compute_resultant_roots(nt_traits,
                                       alg_s, alg_r, alg_t,
                                       alg_v, alg_u, alg_w,
                                       deg1,
                                       arc.m_extra_data->b,
                                       arc.m_extra_data->a,
                                       arc.m_extra_data->c,
                                       ys);
        CGAL_assertion(n_ys <= 2);
      }
    }
    else {
      // Compute the x-coordinates of the intersection points.
      n_xs = compute_resultant_roots(nt_traits,
                                     this->m_r, this->m_s, this->m_t,
                                     this->m_u, this->m_v, this->m_w,
                                     deg1,
                                     arc.m_r, arc.m_s, arc.m_t,
                                     arc.m_u, arc.m_v, arc.m_w,
                                     deg2,
                                     xs);
      CGAL_assertion(n_xs <= 4);

      // Compute the y-coordinates of the intersection points.
      n_ys = compute_resultant_roots(nt_traits,
                                     this->m_s, this->m_r, this->m_t,
                                     this->m_v, this->m_u, this->m_w,
                                     deg1,
                                     arc.m_s, arc.m_r, arc.m_t,
                                     arc.m_v, arc.m_u, arc.m_w,
                                     deg2,
                                     ys);
      CGAL_assertion(n_ys <= 4);
    }

    // Pair the coordinates of the intersection points. As the vectors of
    // x and y-coordinates are sorted in ascending order, we output the
    // intersection points in lexicographically ascending order.
    unsigned int  mult;
    int i, j;

    if (arc.is_special_segment()) {
      if ((n_xs == 0) || (n_ys == 0)) return;

      if ((n_xs == 1) && (n_ys == 1)) {
        // Single intersection.
        Conic_point_2 ip (xs[0], ys[0]);

        ip.set_generating_conic(m_id);
        ip.set_generating_conic (arc.m_id);

        // In case the other curve is of degree 2, this is a tangency point.
        mult = ((deg1 == 1) || is_special_segment()) ? 1 : 2;
        inter_list.push_back(Intersection_point(ip, mult));
      }
      else if ((n_xs == 1) && (n_ys == 2)) {
        Conic_point_2 ip1(xs[0], ys[0]);

        ip1.set_generating_conic(m_id);
        ip1.set_generating_conic(arc.m_id);

        inter_list.push_back(Intersection_point(ip1, 1));

        Conic_point_2 ip2 (xs[0], ys[1]);

        ip2.set_generating_conic(m_id);
        ip2.set_generating_conic(arc.m_id);

        inter_list.push_back(Intersection_point(ip2, 1));
      }
      else if ((n_xs == 2) && (n_ys == 1)) {
        Conic_point_2 ip1 (xs[0], ys[0]);

        ip1.set_generating_conic(m_id);
        ip1.set_generating_conic (arc.m_id);

        inter_list.push_back(Intersection_point(ip1, 1));

        Conic_point_2 ip2 (xs[1], ys[0]);

        ip2.set_generating_conic(m_id);
        ip2.set_generating_conic(arc.m_id);

        inter_list.push_back(Intersection_point(ip2, 1));

      }
      else {
        CGAL_assertion((n_xs == 2) && (n_ys == 2));

        // The x-coordinates and the y-coordinates are given in ascending
        // order. If the slope of the segment is positive, we pair the
        // coordinates as is - otherwise, we swap the pairs.
        int ind_first_y = 0, ind_second_y = 1;

        if (CGAL::sign(arc.m_extra_data->b) == CGAL::sign(arc.m_extra_data->a)) {
          ind_first_y = 1;
          ind_second_y = 0;
        }

        Conic_point_2 ip1(xs[0], ys[ind_first_y]);

        ip1.set_generating_conic(m_id);
        ip1.set_generating_conic(arc.m_id);

        inter_list.push_back(Intersection_point(ip1, 1));

        Conic_point_2 ip2(xs[1], ys[ind_second_y]);

        ip2.set_generating_conic(m_id);
        ip2.set_generating_conic(arc.m_id);

        inter_list.push_back(Intersection_point(ip2, 1));
      }

      return;
    }

    for (i = 0; i < n_xs; ++i) {
      for (j = 0; j < n_ys; ++j) {
        if (is_on_supporting_conic (xs[i], ys[j]) &&
            arc.is_on_supporting_conic(xs[i], ys[j]))
        {
          // Create the intersection point and set its generating conics.
          Conic_point_2 ip(xs[i], ys[j]);

          ip.set_generating_conic(m_id);
          ip.set_generating_conic(arc.m_id);

          // Compute the multiplicity of the intersection point.
          if (deg1 == 1 && deg2 == 1) mult = 1;
          else mult = multiplicity_of_intersection_point(arc, ip);

          // Insert the intersection point to the output list.
          inter_list.push_back(Intersection_point(ip, mult));
        }
      }
    }
  }

  /*! Compute the multiplicity of an intersection point.
   * \param arc The arc to intersect with.
   * \param p The intersection point.
   * \return The multiplicity of the intersection point.
   */
  unsigned int multiplicity_of_intersection_point(const Self& arc,
                                                  const Point_2& p) const {
    CGAL_assertion(! is_special_segment() || ! arc.is_special_segment());

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic slope1_numer, slope1_denom;
    Algebraic slope2_numer, slope2_denom;

    derive_by_x_at(p, 1, slope1_numer, slope1_denom);
    arc.derive_by_x_at(p, 1, slope2_numer, slope2_denom);

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
      arc.derive_by_x_at(p, 2, slope2_numer, slope2_denom);
    }
    else {
      // Both curves have a vertical slope at p.
      // Compare their second-order derivative by y:
      derive_by_y_at(p, 2, slope1_numer, slope1_denom);
      arc.derive_by_y_at(p, 2, slope2_numer, slope2_denom);
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
                         const _Conic_x_monotone_arc_2<Conic_arc_2>& arc)
{
  // Output the supporting conic curve.
  os << "{" << CGAL::to_double(arc.r()) << "*x^2 + "
     << CGAL::to_double(arc.s()) << "*y^2 + "
     << CGAL::to_double(arc.t()) << "*xy + "
     << CGAL::to_double(arc.u()) << "*x + "
     << CGAL::to_double(arc.v()) << "*y + "
     << CGAL::to_double(arc.w()) << "}";

  // Output the endpoints.
  os << " : (" << CGAL::to_double(arc.source().x()) << ","
     << CGAL::to_double(arc.source().y()) << ") ";

  if (arc.orientation() == CLOCKWISE) os << "--cw-->";
  else if (arc.orientation() == COUNTERCLOCKWISE) os << "--ccw-->";
  else os << "--l-->";

  os << " (" << CGAL::to_double(arc.target().x()) << ","
     << CGAL::to_double(arc.target().y()) << ")";

  return (os);
}

} //namespace CGAL

#endif
