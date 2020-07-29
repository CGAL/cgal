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

template <typename Conic_arc_>
class _Conic_x_monotone_arc_2 : private Conic_arc_ {
public:

  typedef Conic_arc_                              Conic_arc_2;
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

  using Conic_arc_2::_sign_of_extra_data;
  using Conic_arc_2::_is_between_endpoints;
  using Conic_arc_2::_is_strictly_between_endpoints;
  using Conic_arc_2::_conic_get_y_coordinates;
  /*!
   * \struct Less functor for Conic_pair.
   */
  struct Less_conic_pair {
    bool operator()(const Conic_pair& cp1, const Conic_pair& cp2) const
    {
      // Compare the pairs of IDs lexicographically.
      return ((cp1.first < cp2.first) ||
              ((cp1.first == cp2.first) && (cp1.second < cp2.second)));
    }
  };

  typedef std::map<Conic_pair,
                   Intersection_list,
                   Less_conic_pair>               Intersection_map;
  typedef typename Intersection_map::value_type   Intersection_map_entry;
  typedef typename Intersection_map::iterator     Intersection_map_iterator;

protected:

  typedef Conic_arc_2                             Base;

  typedef typename Conic_arc_2::Integer           Integer;
  typedef typename Conic_arc_2::Nt_traits         Nt_traits;
  typedef typename Conic_arc_2::Rat_kernel        Rat_kernel;

  // Bit masks for the _info field (the two least significant bits are already
  // used by the base class).
  enum {
    IS_VERTICAL_SEGMENT = 4,
    IS_DIRECTED_RIGHT = 8,
    DEGREE_1 = 16,
    DEGREE_2 = 32,
    DEGREE_MASK = 16 + 32,
    PLUS_SQRT_DISC_ROOT = 64,
    FACING_UP = 128,
    FACING_DOWN = 256,
    FACING_MASK = 128 + 256,
    IS_SPECIAL_SEGMENT = 512
  };

  Algebraic alg_r;      // The coefficients of the supporting conic curve:
  Algebraic alg_s;      //
  Algebraic alg_t;      //   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0 ,
  Algebraic alg_u;      //
  Algebraic alg_v;      // converted to algebraic numbers.
  Algebraic alg_w;      //

  Conic_id _id;         // The ID number of the supporting conic curve.

public:

  /// \name Constrcution methods.
  //@{

  /*!
   * Default constructor.
   */
  _Conic_x_monotone_arc_2 () :
    Base(),
    _id()
  {}

  /*!
   * Copy constructor.
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
    _id(arc._id)
  {}

  /*! Construct an x-monotone arc from a conic arc.
   * \param arc The given (base) arc.
   * \pre The given arc is x-monotone.
   */
  _Conic_x_monotone_arc_2(const Base& arc) :
    Base(arc),
    _id()
  {
    CGAL_precondition(arc.is_valid() && arc.is_x_monotone());
    _set ();
  }

  /*! Construct an x-monotone arc from a conic arc.
   * \param arc The given (base) arc.
   * \param id The ID of the base arc.
   */
  _Conic_x_monotone_arc_2(const Base& arc, const Conic_id& id) :
    Base(arc),
    _id(id)
  {
    CGAL_precondition(arc.is_valid() && id.is_valid());
    _set ();
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
    _id(id)
  {
    CGAL_precondition(arc.is_valid() && id.is_valid());

    // Set the two endpoints.
    this->_source = source;
    this->_target = target;

    _set();
  }

  /*!
   * Construct a special segment connecting to given endpoints (for the usage
   * of the landmarks point-location strategy).
   * \param source The source point.
   * \param target The target point.
   */
  _Conic_x_monotone_arc_2(const Point_2& source, const Point_2& target) :
    Base()
  {
    // Set the basic properties and clear the _info bits.
    this->_source = source;
    this->_target = target;
    this->_orient = COLLINEAR;
    this->_info = 0;

    // Check if the arc is directed right (the target is lexicographically
    // greater than the source point), or to the left.
    Alg_kernel ker;
    Comparison_result dir_res =
      ker.compare_xy_2_object()(this->_source, this->_target);

    CGAL_precondition (dir_res != EQUAL);
    // Invalid arc:
    if (dir_res == EQUAL) return;

    this->_info = (Conic_arc_2::IS_VALID | DEGREE_1);
    if (dir_res == SMALLER)
      this->_info = (this->_info | IS_DIRECTED_RIGHT);

    // Compose the equation of the underlying line.
    const Algebraic x1 = source.x(), y1 = source.y();
    const Algebraic x2 = target.x(), y2 = target.y();

    // The supporting line is A*x + B*y + C = 0, where:
    //
    //  A = y2 - y1,    B = x1 - x2,    C = x2*y1 - x1*y2
    //
    // We use the extra data field to store the equation of this line.
    this->_extra_data_P = new typename Base::Extra_data;
    this->_extra_data_P->a = y2 - y1;
    this->_extra_data_P->b = x1 - x2;
    this->_extra_data_P->c = x2*y1 - x1*y2;
    this->_extra_data_P->side = ZERO;

    // Check if the segment is vertical.
    if (CGAL::sign (this->_extra_data_P->b) == ZERO)
      this->_info = (this->_info | IS_VERTICAL_SEGMENT);

    // Mark that this is a special segment.
    this->_info = (this->_info | IS_SPECIAL_SEGMENT);

    return;
  }

  /*!
   * Construct a special segment of a given line connecting to given
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
    CGAL_precondition(CGAL::sign(a * source.x() +
                                 b * source.y() + c) == CGAL::ZERO);

    CGAL_precondition(CGAL::sign(a * target.x() +
                                 b * target.y() + c) == CGAL::ZERO);

    // Set the basic properties and clear the _info bits.
    this->_source = source;
    this->_target = target;
    this->_orient = COLLINEAR;
    this->_info = 0;

    // Check if the arc is directed right (the target is lexicographically
    // greater than the source point), or to the left.
    Alg_kernel ker;
    Comparison_result res =
      ker.compare_x_2_object()(this->_source, this->_target);

    this->_info = (Conic_arc_2::IS_VALID | DEGREE_1);
    if (res == EQUAL) {
      // Mark that the segment is vertical.
      this->_info = (this->_info | IS_VERTICAL_SEGMENT);

      // Compare the endpoints lexicographically.
      res = ker.compare_y_2_object()(this->_source, this->_target);

      CGAL_precondition (res != EQUAL);
      if (res == EQUAL) {
        // Invalid arc:
        this->_info = 0;
        return;
      }
    }

    if (res == SMALLER)
      this->_info = (this->_info | IS_DIRECTED_RIGHT);

    // Store the coefficients of the line.
    this->_extra_data_P = new typename Base::Extra_data;
    this->_extra_data_P->a = a;
    this->_extra_data_P->b = b;
    this->_extra_data_P->c = c;
    this->_extra_data_P->side = ZERO;

    // Mark that this is a special segment.
    this->_info = (this->_info | IS_SPECIAL_SEGMENT);

    return;
  }

  /*!
   * Assignment operator.
   * \param arc The copied arc.
   */
  const Self& operator=(const Self& arc)
  {
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

    _id = arc._id;

    return (*this);
  }
  //@}

  /// \name Accessing the arc properties.
  //@{

  /*!
   * Get the coefficients of the underlying conic.
   */
  const Integer& r() const { return (this->_r); }
  const Integer& s() const { return (this->_s); }
  const Integer& t() const { return (this->_t); }
  const Integer& u() const { return (this->_u); }
  const Integer& v() const { return (this->_v); }
  const Integer& w() const { return (this->_w); }

  /*!
   * Get the arc's source.
   * \return The source point.
   */
  const Conic_point_2& source() const { return (this->_source); }

  /*! Get the arc's target.
   * \return The target point.
   */
  const Conic_point_2& target() const { return (this->_target); }

  /*! Get the orientation of the arc.
   * \return The orientation.
   */
  Orientation orientation() const { return (this->_orient); }

  /*! Get the left endpoint of the arc.
   */
  const Conic_point_2& left () const
  {
    if ((this->_info & IS_DIRECTED_RIGHT) != 0) return (this->_source);
    else return (this->_target);
  }

  /*!
   * Get the right endpoint of the arc.
   */
  const Conic_point_2& right () const
  {
    if ((this->_info & IS_DIRECTED_RIGHT) != 0) return (this->_target);
    else return (this->_source);
  }

  /*!
   * Return true iff the conic arc is directed right iexicographically.
   */
  bool is_directed_right() const
  { return ((this->_info & IS_DIRECTED_RIGHT) != 0); }

  /*!
   * Get a bounding box for the conic arc.
   * \return The bounding box.
   */
  Bbox_2 bbox() const { return (Base::bbox()); }
  //@}

  /// \name Predicates.
  //@{

  /*!
   * Check if the conic arc is a vertical segment.
   */
  bool is_vertical () const
  {
    return ((this->_info & IS_VERTICAL_SEGMENT) != 0);
  }

  /*!
   * Check whether the given point lies on the arc.
   * \param p The qury point.
   * \param (true) if p lies on the arc; (false) otherwise.
   */
  bool contains_point(const Conic_point_2& p) const
  {
    // First check if p lies on the supporting conic. We first check whether
    // it is one of p's generating conic curves.
    bool       p_on_conic = false;

    if (p.is_generating_conic(_id)) {
      p_on_conic = true;
    }
    else {
      // Check whether p satisfies the supporting conic equation.
      p_on_conic = _is_on_supporting_conic(p.x(), p.y());

      if (p_on_conic) {
        // As p lies on the supporting conic of our arc, add its ID to
        // the list of generating conics for p.
        Conic_point_2&  p_non_const = const_cast<Conic_point_2&> (p);
        p_non_const.set_generating_conic (_id);
      }
    }

    if (! p_on_conic)
      return (false);

    // Check if p is between the endpoints of the arc.
    return (_is_between_endpoints (p));
  }
  //@}

  /// \name Constructing points on the arc.
  //@{

  /*!
   * Compute a point on the arc with the same x-coordiante as the given point.
   * \param p The given point.
   * \pre The arc is not vertical and p is in the x-range of the arc.
   * \return A point on the arc with the same x-coordiante as p.
   */
  Point_2 point_at_x (const Point_2& p) const
  {
    // Make sure that p is in the x-range of the arc.
    CGAL_precondition ((this->_info & IS_VERTICAL_SEGMENT) == 0);

    CGAL_precondition_code (
      Alg_kernel   ker;
    );

    CGAL_precondition(ker.compare_x_2_object() (p, left()) != SMALLER &&
                      ker.compare_x_2_object() (p, right()) != LARGER);

    if (_is_special_segment()) {
      // In case of a special segment, the equation of the supported line
      // (a*x + b*y + c) = 0 is stored with the extra data field, and we
      // simply have:
      Algebraic _y = -(this->_extra_data_P->a*p.x() +
                       this->_extra_data_P->c) /
        this->_extra_data_P->b;

      // Return the computed point.
      return (Point_2 (p.x(), _y));
    }

    // Compute the y-coordinate according to the degree of the supporting
    // conic curve.
    Nt_traits nt_traits;
    Algebraic y;

    if ((this->_info & DEGREE_MASK) == DEGREE_1) {
      // In case of a linear curve, the y-coordinate is a simple linear
      // expression of x(p) (note that v is not 0 as the arc is not vertical):
      //   y = -(u*x(p) + w) / v
      y = -(alg_u*p.x() + alg_w) / alg_v;
    }
    else if (this->_orient == COLLINEAR) {
      CGAL_assertion (this->_extra_data_P != nullptr);

      // In this case the equation of the supporting line is given by the
      // extra data structure.
      y = -(this->_extra_data_P->a * p.x() +
            this->_extra_data_P->c) / this->_extra_data_P->b;
    }
    else {
      CGAL_assertion((this->_info & DEGREE_MASK) == DEGREE_2);

      // In this case the y-coordinate is one of solutions to the quadratic
      // equation:
      //  s*y^2 + (t*x(p) + v)*y + (r*x(p)^2 + u*x(p) + w) = 0
      Algebraic A = alg_s;
      Algebraic B = alg_t*p.x() + alg_v;
      Algebraic C = (alg_r*p.x() + alg_u)*p.x() + alg_w;

      if (CGAL::sign(this->_s) == ZERO) {
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
        if ((this->_info & PLUS_SQRT_DISC_ROOT) != 0) {
          y = (nt_traits.sqrt (disc) - B) / (2*A);
        }
        else {
          y = -(B + nt_traits.sqrt (disc)) / (2*A);
        }
      }
    }

    // Return the computed point.
    return (Point_2 (p.x(), y));
  }

  /*!
   * Get a polyline approximating the conic arc.
   * \param n The maximal number of sample points.
   * \param oi An output iterator, whose value-type is pair<double,double>
   *           (representing an approximated point).
   *           In case the arc is a line segment, there are 2 output points,
   *           otherwise the arc is approximated by the polyline defined by
   *           (p_0, p_1, ..., p_n), where p_0 and p_n are the left and right
   *           endpoints of the arc, respectively.
   */
  template <class OutputIterator>
  OutputIterator polyline_approximation (size_t n,
                                         OutputIterator oi) const
  {
    CGAL_precondition (n != 0);

    const double x_left = CGAL::to_double (left().x());
    const double y_left = CGAL::to_double (left().y());
    const double x_right = CGAL::to_double (right().x());
    const double y_right = CGAL::to_double (right().y());

    if (this->_orient == COLLINEAR) {
      // In case of a line segment, return the two endpoints.
      *oi++ = std::pair<double, double> (x_left, y_left);
      *oi++ = std::pair<double, double> (x_right, y_right);
      return oi;
    }

    // Otherwise, sample (n - 1) equally-spaced points in between.
    const double app_r = CGAL::to_double (this->_r);
    const double app_s = CGAL::to_double (this->_s);
    const double app_t = CGAL::to_double (this->_t);
    const double app_u = CGAL::to_double (this->_u);
    const double app_v = CGAL::to_double (this->_v);
    const double app_w = CGAL::to_double (this->_w);
    const double x_jump = (x_right - x_left) / n;
    double       x, y;
    const bool   A_is_zero = (CGAL::sign(this->_s) == ZERO);
    double       A = app_s, B, C;
    double       disc;
    size_t       i;

    *oi = std::pair<double, double>(x_left, y_left);   // The left point.
    ++oi;
    for (i = 1; i < n; i++) {
      x = x_left + x_jump*i;

      // Solve the quadratic equation: A*x^2 + B*x + C = 0:
      B = app_t*x + app_v;
      C = (app_r*x + app_u)*x + app_w;

      if (A_is_zero) {
        y = -C / B;
      }
      else {
        disc = B*B - 4*A*C;

        if (disc < 0) disc = 0;

        // We take either the root involving -sqrt(disc) or +sqrt(disc)
        // based on the information flags.
        if ((this->_info & PLUS_SQRT_DISC_ROOT) != 0) {
          y = (std::sqrt(disc) - B) / (2*A);
        }
        else {
          y = -(B + std::sqrt (disc)) / (2*A);
        }
      }

      *oi++ = std::pair<double, double> (x, y);
    }
    *oi++ = std::pair<double, double> (x_right, y_right);   // The right point.

    return oi;
  }

  /*! Compare to arcs immediately to the right of their intersection point.
   * \param arc The compared arc.
   * \param p The reference intersection point.
   * \return The relative position of the arcs to the right of p.
   * \pre Both arcs we compare are not vertical segments.
   */
  Comparison_result compare_to_right(const Self& arc,
                                     const Conic_point_2& p) const
  {
    CGAL_precondition((this->_info & IS_VERTICAL_SEGMENT) == 0 &&
                      (arc._info & IS_VERTICAL_SEGMENT) == 0);

    // In case one arc is facing upwards and another facing downwards, it is
    // clear that the one facing upward is above the one facing downwards.
    if (_has_same_supporting_conic (arc)) {
      if ((this->_info & FACING_UP) != 0 && (arc._info & FACING_DOWN) != 0)
        return LARGER;
      else if ((this->_info & FACING_DOWN)!= 0 && (arc._info & FACING_UP) != 0)
        return SMALLER;

      // In this case the two arcs overlap.
      CGAL_assertion((this->_info & FACING_MASK) == (arc._info & FACING_MASK));

      return EQUAL;
    }

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic slope1_numer, slope1_denom;
    Algebraic slope2_numer, slope2_denom;

    _derive_by_x_at (p, 1, slope1_numer, slope1_denom);
    arc._derive_by_x_at (p, 1, slope2_numer, slope2_denom);

    // Check if any of the slopes is vertical.
    const bool is_vertical_slope1 = (CGAL::sign(slope1_denom) == ZERO);
    const bool is_vertical_slope2 = (CGAL::sign(slope2_denom) == ZERO);

    if (!is_vertical_slope1 && !is_vertical_slope2) {
      // The two derivatives at p are well-defined: use them to determine
      // which arc is above the other (the one with a larger slope is below).
      Comparison_result slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the second-order derivative.
      _derive_by_x_at(p, 2, slope1_numer, slope1_denom);
      arc._derive_by_x_at(p, 2, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the third-order derivative.
      _derive_by_x_at(p, 3, slope1_numer, slope1_denom);
      arc._derive_by_x_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      // \todo Handle higher-order derivatives:
      CGAL_assertion (slope_res != EQUAL);

      return (slope_res);
    }
    else if (!is_vertical_slope2) {
      // The first arc has a vertical slope at p: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion ((this->_info & FACING_MASK) != 0);

      if ((this->_info & FACING_UP) != 0) return (LARGER);
      return SMALLER;
    }
    else if (!is_vertical_slope1) {
      // The second arc has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion ((arc._info & FACING_MASK) != 0);

      if ((arc._info & FACING_UP) != 0) return (SMALLER);
      return LARGER;
    }

    // The two arcs have vertical slopes at p_int:
    // First check whether one is facing up and one down. In this case the
    // comparison result is trivial.
    if ((this->_info & FACING_UP) != 0 && (arc._info & FACING_DOWN) != 0)
      return (LARGER);
    else if ((this->_info & FACING_DOWN)!= 0 && (arc._info & FACING_UP)!= 0)
      return SMALLER;

    // Compute the second-order derivative by y and act according to it.
    _derive_by_y_at (p, 2, slope1_numer, slope1_denom);
    arc._derive_by_y_at (p, 2, slope2_numer, slope2_denom);

    Comparison_result slope_res =
      CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

    // If necessary, use the third-order derivative by y.
    if (slope_res == EQUAL) {
      // \todo Check this!
      _derive_by_y_at(p, 3, slope1_numer, slope1_denom);
      arc._derive_by_y_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);
    }

    // \todo Handle higher-order derivatives:
    CGAL_assertion(slope_res != EQUAL);

    if ((this->_info & FACING_UP) != 0 && (arc._info & FACING_UP) != 0) {
      // Both are facing up.
      return ((slope_res == LARGER) ? SMALLER : LARGER);
    }
    // Both are facing down.
    return (slope_res);
  }

  /*!
   * Compare to arcs immediately to the leftt of their intersection point.
   * \param arc The compared arc.
   * \param p The reference intersection point.
   * \return The relative position of the arcs to the left of p.
   * \pre Both arcs we compare are not vertical segments.
   */
  Comparison_result compare_to_left(const Self& arc,
                                    const Conic_point_2& p) const
  {
    CGAL_precondition((this->_info & IS_VERTICAL_SEGMENT) == 0 &&
                      (arc._info & IS_VERTICAL_SEGMENT) == 0);

    // In case one arc is facing upwards and another facing downwards, it is
    // clear that the one facing upward is above the one facing downwards.
    if (_has_same_supporting_conic (arc)) {
      if ((this->_info & FACING_UP) != 0 && (arc._info & FACING_DOWN) != 0)
        return LARGER;
      else if ((this->_info & FACING_DOWN)!= 0 && (arc._info & FACING_UP)!= 0)
        return SMALLER;

      // In this case the two arcs overlap.
      CGAL_assertion((this->_info & FACING_MASK) == (arc._info & FACING_MASK));

      return EQUAL;
    }

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic slope1_numer, slope1_denom;
    Algebraic slope2_numer, slope2_denom;

    _derive_by_x_at(p, 1, slope1_numer, slope1_denom);
    arc._derive_by_x_at(p, 1, slope2_numer, slope2_denom);

    // Check if any of the slopes is vertical.
    const bool is_vertical_slope1 = (CGAL::sign (slope1_denom) == ZERO);
    const bool is_vertical_slope2 = (CGAL::sign (slope2_denom) == ZERO);

    if (!is_vertical_slope1 && !is_vertical_slope2) {
      // The two derivatives at p are well-defined: use them to determine
      // which arc is above the other (the one with a larger slope is below).
      Comparison_result  slope_res = CGAL::compare(slope2_numer*slope1_denom,
                                                   slope1_numer*slope2_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the second-order derivative.
      _derive_by_x_at (p, 2, slope1_numer, slope1_denom);
      arc._derive_by_x_at (p, 2, slope2_numer, slope2_denom);

      slope_res = CGAL::compare (slope1_numer*slope2_denom,
                                 slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the third-order derivative.
      _derive_by_x_at(p, 3, slope1_numer, slope1_denom);
      arc._derive_by_x_at(p, 3, slope2_numer, slope2_denom);

      slope_res = CGAL::compare(slope2_numer*slope1_denom,
                                slope1_numer*slope2_denom);

      // \todo Handle higher-order derivatives:
      CGAL_assertion (slope_res != EQUAL);

      return (slope_res);
    }
    else if (!is_vertical_slope2) {
      // The first arc has a vertical slope at p: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion ((this->_info & FACING_MASK) != 0);

      if ((this->_info & FACING_UP) != 0) return (LARGER);
      return SMALLER;
    }
    else if (!is_vertical_slope1) {
      // The second arc has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion ((arc._info & FACING_MASK) != 0);

      if ((arc._info & FACING_UP) != 0) return (SMALLER);
      return LARGER;
    }

    // The two arcs have vertical slopes at p_int:
    // First check whether one is facing up and one down. In this case the
    // comparison result is trivial.
    if ((this->_info & FACING_UP) != 0 && (arc._info & FACING_DOWN) != 0)
      return LARGER;
    else if ((this->_info & FACING_DOWN)!= 0 && (arc._info & FACING_UP)!= 0)
      return SMALLER;

    // Compute the second-order derivative by y and act according to it.
    _derive_by_y_at(p, 2, slope1_numer, slope1_denom);
    arc._derive_by_y_at(p, 2, slope2_numer, slope2_denom);

    Comparison_result slope_res =
      CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);

    // If necessary, use the third-order derivative by y.
    if (slope_res == EQUAL) {
      // \todo Check this!
      _derive_by_y_at(p, 3, slope1_numer, slope1_denom);
      arc._derive_by_y_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);
    }

    // \todo Handle higher-order derivatives:
    CGAL_assertion(slope_res != EQUAL);

    if ((this->_info & FACING_UP) != 0 && (arc._info & FACING_UP) != 0) {
      // Both are facing up.
      return ((slope_res == LARGER) ? SMALLER : LARGER);
    }
    // Both are facing down.
    return (slope_res);
  }

  /*!
   * Compute the intersections with the given arc.
   * \param arc The given intersecting arc.
   * \param inter_map Maps conic pairs to lists of their intersection points.
   * \param oi The output iterator.
   * \return The past-the-end iterator.
   */
  template <typename OutputIterator>
  OutputIterator intersect(const Self& arc,
                           Intersection_map& inter_map,
                           OutputIterator oi) const
  {
    typedef boost::variant<Intersection_point, Self>  Intersection_result;

    if (_has_same_supporting_conic(arc)) {
      // Check for overlaps between the two arcs.
      Self overlap;

      if (_compute_overlap(arc, overlap)) {
        // There can be just a single overlap between two x-monotone arcs:
        *oi++ = Intersection_result(overlap);
        return oi;
      }

      // In case there is not overlap and the supporting conics are the same,
      // there cannot be any intersection points, unless the two arcs share
      // an end point.
      // Note that in this case we do not define the multiplicity of the
      // intersection points we report.
      Alg_kernel  ker;

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

    if (_id.is_valid() && arc._id.is_valid()) {
      if (_id < arc._id) conic_pair = Conic_pair (_id, arc._id);
      else conic_pair = Conic_pair (arc._id, _id);
      map_iter = inter_map.find (conic_pair);
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
      _intersect_supporting_conics(arc, inter_list);

      if (! invalid_ids) inter_map[conic_pair] = inter_list;
    }
    else {
      // Obtain the precomputed intersection points from the map.
      inter_list = (*map_iter).second;
    }

    // Go over the list of intersection points and report those that lie on
    // both x-monotone arcs.
    for (auto iter = inter_list.begin(); iter != inter_list.end(); ++iter) {
      if (_is_between_endpoints((*iter).first) &&
          arc._is_between_endpoints((*iter).first))
      {
        *oi++ = Intersection_result(*iter);
      }
    }

    return oi;
  }
  //@}

  /// \name Constructing x-monotone arcs.
  //@{

  /*!
   * Split the arc into two at a given split point.
   * \param p The split point.
   * \param c1 Output: The first resulting arc, lying to the left of p.
   * \param c2 Output: The first resulting arc, lying to the right of p.
   * \pre p lies in the interior of the arc (not one of its endpoints).
   */
  void split(const Conic_point_2& p, Self& c1, Self& c2) const
  {
    // Make sure that p lies on the interior of the arc.
    CGAL_precondition_code(Alg_kernel ker);

    CGAL_precondition (this->contains_point (p) &&
                       ! ker.equal_2_object() (p, this->_source) &&
                       ! ker.equal_2_object() (p, this->_target));

    // Make copies of the current arc.
    c1 = *this;
    c2 = *this;

    // Assign the endpoints of the arc.
    if ((this->_info & IS_DIRECTED_RIGHT) != 0)
    {
      // The arc is directed from left to right, so p becomes c1's target
      // and c2's source.
      c1._target = p;
      c2._source = p;

      if (! p.is_generating_conic (_id)) {
        c1._target.set_generating_conic (_id);
        c2._source.set_generating_conic (_id);
      }
    }
    else
    {
      // The arc is directed from right to left, so p becomes c2's target
      // and c1's source.
      c1._source = p;
      c2._target = p;

      if (! p.is_generating_conic (_id)) {
        c1._source.set_generating_conic (_id);
        c2._target.set_generating_conic (_id);
      }
    }

    return;
  }

  /*!
   * Flip the arc.
   * \return An arc with swapped source and target and a reverse orienation.
   */
  Self flip() const
  {
    // Make a copy of the current arc.
    Self arc = *this;

    // Reverse the orientation.
    if (this->_orient == CLOCKWISE) arc._orient = COUNTERCLOCKWISE;
    else if (this->_orient == COUNTERCLOCKWISE) arc._orient = CLOCKWISE;

    // Swap the source and the target.
    arc._source = this->_target;
    arc._target = this->_source;

    // Change the direction bit among the information flags.
    arc._info = (this->_info ^ IS_DIRECTED_RIGHT);

    return arc;
  }

  /*!
   * Trim the arc given its new endpoints.
   * \param ps The new source point.
   * \param pt The new target point.
   * \return The new trimmed arc.
   * \pre Both ps and pt lies on the arc and must conform with the current
   *      direction of the arc.
   */
  Self trim(const Conic_point_2& ps, const Conic_point_2& pt) const
  {
    // Make sure that both ps and pt lie on the arc.
    CGAL_precondition(this->contains_point (ps) &&
                      this->contains_point (pt));

    // Make sure that the endpoints conform with the direction of the arc.
    Self         arc = *this;
    Alg_kernel   ker;

    if (! ((((this->_info & IS_DIRECTED_RIGHT) != 0) &&
            ker.compare_xy_2_object() (ps, pt) == SMALLER) ||
           (((this->_info & IS_DIRECTED_RIGHT) == 0) &&
            ker.compare_xy_2_object() (ps, pt) == LARGER)))
    {
      // We are allowed to change the direction only in case of a segment.
      CGAL_assertion (this->_orient == COLLINEAR);
      arc._info = (this->_info ^ IS_DIRECTED_RIGHT);
    }

    // Make a copy of the current arc and assign its endpoints.
    if (! ker.equal_2_object() (ps, this->_source)) {
      arc._source = ps;

      if (! ps.is_generating_conic (_id))
        arc._source.set_generating_conic (_id);
    }

    if (! ker.equal_2_object() (pt, this->_target))
    {
      arc._target = pt;

      if (! pt.is_generating_conic (_id))
        arc._target.set_generating_conic (_id);
    }

    return (arc);
  }

  /*!
   * Check whether the two arcs are equal (have the same graph).
   * \param arc The compared arc.
   * \return (true) if the two arcs have the same graph; (false) otherwise.
   */
  bool equals (const Self& arc) const
  {
    // The two arc must have the same supporting conic curves.
    if (! _has_same_supporting_conic (arc))
      return false;

    // Check that the arc endpoints are the same.
    Alg_kernel ker;

    if (this->_orient == COLLINEAR) {
      CGAL_assertion(arc._orient == COLLINEAR);
      return((ker.equal_2_object()(this->_source, arc._source) &&
              ker.equal_2_object()(this->_target, arc._target)) ||
              (ker.equal_2_object()(this->_source, arc._target) &&
               ker.equal_2_object()(this->_target, arc._source)));
    }

    if (this->_orient == arc._orient) {
      // Same orientation - the source and target points must be the same.
      return (ker.equal_2_object()(this->_source, arc._source) &&
              ker.equal_2_object()(this->_target, arc._target));
    }
    else {
      // Reverse orientation - the source and target points must be swapped.
      return (ker.equal_2_object()(this->_source, arc._target) &&
              ker.equal_2_object()(this->_target, arc._source));
    }
  }

  /*! Check whether it is possible to merge the arc with the given arc.
   * \param arc The query arc.
   * \return (true) if it is possible to merge the two arcs;
   *         (false) otherwise.
   */
  bool can_merge_with(const Self& arc) const
  {
    // In order to merge the two arcs, they should have the same supporting
    // conic.
    if (! _has_same_supporting_conic(arc)) return false;

    // Check if the left endpoint of one curve is the right endpoint of the
    // other.
    Alg_kernel   ker;

    return (ker.equal_2_object() (right(), arc.left()) ||
            ker.equal_2_object() (left(), arc.right()));
  }

  /*! Merge the current arc with the given arc.
   * \param arc The arc to merge with.
   * \pre The two arcs are mergeable.
   */
  void merge(const Self& arc)
  {
    CGAL_precondition (this->can_merge_with (arc));

    // Check if we should extend the arc to the left or to the right.
    Alg_kernel   ker;

    if (ker.equal_2_object() (right(), arc.left())) {
      // Extend the arc to the right.
      if ((this->_info & IS_DIRECTED_RIGHT) != 0) this->_target = arc.right();
      else this->_source = arc.right();
    }
    else {
      CGAL_precondition (ker.equal_2_object() (left(), arc.right()));

      // Extend the arc to the left.
      if ((this->_info & IS_DIRECTED_RIGHT) != 0)
        this->_source = arc.left();
      else
        this->_target = arc.left();
    }

    return;
  }

  bool is_upper() const
  {
    return ((this->_info & FACING_UP) != 0);
  }

  bool is_lower() const
  {
    return ((this->_info & FACING_DOWN) != 0);
  }
  //@}

private:

  /// \name Auxiliary (private) functions.
  //@{

  /*!
   * Set the properties of the x-monotone conic arc (for the usage of the
   * constructors).
   */
  void _set ()
  {
    // Convert the coefficients of the supporting conic to algebraic numbers.
    Nt_traits        nt_traits;

    alg_r = nt_traits.convert (this->_r);
    alg_s = nt_traits.convert (this->_s);
    alg_t = nt_traits.convert (this->_t);
    alg_u = nt_traits.convert (this->_u);
    alg_v = nt_traits.convert (this->_v);
    alg_w = nt_traits.convert (this->_w);

    // Set the generating conic ID for the source and target points.
    this->_source.set_generating_conic (_id);
    this->_target.set_generating_conic (_id);

    // Clear the _info bits.
    this->_info = Conic_arc_2::IS_VALID;

    // Check if the arc is directed right (the target is lexicographically
    // greater than the source point), or to the left.
    Alg_kernel         ker;
    Comparison_result  dir_res = ker.compare_xy_2_object() (this->_source,
                                                            this->_target);

    CGAL_assertion (dir_res != EQUAL);

    if (dir_res == SMALLER)
      this->_info = (this->_info | IS_DIRECTED_RIGHT);

    // Compute the degree of the underlying conic.
    if (CGAL::sign (this->_r) != ZERO ||
        CGAL::sign (this->_s) != ZERO ||
        CGAL::sign (this->_t) != ZERO)
    {
      this->_info = (this->_info | DEGREE_2);

      if (this->_orient == COLLINEAR)
      {
        this->_info = (this->_info | IS_SPECIAL_SEGMENT);

        if (ker.compare_x_2_object() (this->_source, this->_target) == EQUAL)
        {
          // The arc is a vertical segment:
          this->_info = (this->_info | IS_VERTICAL_SEGMENT);
        }

        return;
      }
    }
    else
    {
      CGAL_assertion (CGAL::sign (this->_u) != ZERO ||
                      CGAL::sign (this->_v) != ZERO);

      if (CGAL::sign (this->_v) == ZERO)
      {

        // The supporting curve is of the form: _u*x + _w = 0
        this->_info = (this->_info | IS_VERTICAL_SEGMENT);
      }

      this->_info = (this->_info | DEGREE_1);

      return;
    }

    if (this->_orient == COLLINEAR)
      return;

    // Compute a midpoint between the source and the target and get the y-value
    // of the arc at its x-coordiante.
    Point_2          p_mid = ker.construct_midpoint_2_object() (this->_source,
                                                                this->_target);
    Algebraic        ys[2];
    CGAL_assertion_code(int              n_ys = )
      _conic_get_y_coordinates (p_mid.x(), ys);

    CGAL_assertion (n_ys != 0);

    // Check which solution lies on the x-monotone arc.
    Point_2          p_arc_mid (p_mid.x(), ys[0]);

    if (_is_strictly_between_endpoints (p_arc_mid))
    {
      // Mark that we should use the -sqrt(disc) root for points on this
      // x-monotone arc.
      this->_info = (this->_info & ~PLUS_SQRT_DISC_ROOT);
    }
    else
    {
      CGAL_assertion (n_ys == 2);
      p_arc_mid = Point_2 (p_mid.x(), ys[1]);

      CGAL_assertion (_is_strictly_between_endpoints (p_arc_mid));

      // Mark that we should use the +sqrt(disc) root for points on this
      // x-monotone arc.
      this->_info = (this->_info | PLUS_SQRT_DISC_ROOT);
    }

    // Check whether the conic is facing up or facing down:
    // Check whether the arc (which is x-monotone of degree 2) lies above or
    // below the segement that contects its two end-points (x1,y1) and (x2,y2).
    // To do that, we find the y coordinate of a point on the arc whose x
    // coordinate is (x1+x2)/2 and compare it to (y1+y2)/2.
    Comparison_result res = ker.compare_y_2_object() (p_arc_mid, p_mid);

    if (res == LARGER)
    {
      // The arc is above the connecting segment, so it is facing upwards.
      this->_info = (this->_info | FACING_UP);
    }
    else if (res == SMALLER)
    {
      // The arc is below the connecting segment, so it is facing downwards.
      this->_info = (this->_info | FACING_DOWN);
    }

    return;
  }

  /*!
   * Check if the arc is a special segment connecting two algebraic endpoints
   * (and has no undelying integer conic coefficients).
   */
  bool _is_special_segment () const
  {
    return ((this->_info & IS_SPECIAL_SEGMENT) != 0);
  }

  /*!
   * Check whether the given point lies on the supporting conic of the arc.
   * \param px The x-coordinate of query point.
   * \param py The y-coordinate of query point.
   * \return (true) if p lies on the supporting conic; (false) otherwise.
   */
  bool _is_on_supporting_conic (const Algebraic& px,
                                const Algebraic& py) const
  {
    CGAL::Sign       _sign;

    if (! _is_special_segment())
    {
      // Check whether p satisfies the conic equation.
      // The point must satisfy: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0.
      _sign = CGAL::sign ((alg_r*px + alg_t*py + alg_u) * px +
                          (alg_s*py + alg_v) * py +
                          alg_w);
    }
    else
    {
      // Check whether p satisfies the equation of the line stored with the
      // extra data.
      _sign = _sign_of_extra_data (px, py);
    }

    return (_sign == ZERO);
  }

  /*!
   * Check whether the two arcs have the same supporting conic.
   * \param arc The compared arc.
   * \return (true) if the two supporting conics are the same.
   */
  bool _has_same_supporting_conic (const Self& arc) const
  {
    // Check if the two arcs originate from the same conic:
    if (_id == arc._id && _id.is_valid() && arc._id.is_valid())
      return (true);

    // In case both arcs are collinear, check if they have the same
    // supporting lines.
    if (this->_orient == COLLINEAR && arc._orient == COLLINEAR)
    {
      // Construct the two supporting lines and compare them.
      Alg_kernel                             ker;
      typename Alg_kernel::Construct_line_2  construct_line =
                                                 ker.construct_line_2_object();
      typename Alg_kernel::Line_2          l1 = construct_line (this->_source,
                                                                this->_target);
      typename Alg_kernel::Line_2          l2 = construct_line (arc._source,
                                                                arc._target);
      typename Alg_kernel::Equal_2         equal = ker.equal_2_object();

      if (equal (l1, l2))
        return (true);

      // Try to compare l1 with the opposite of l2.
      l2 = construct_line (arc._target, arc._source);

      return (equal (l1, l2));
    }
    else if (this->_orient == COLLINEAR || arc._orient == COLLINEAR)
    {
      // Only one arc is collinear, so the supporting curves cannot be the
      // same:
      return (false);
    }

    // Check whether the coefficients of the two supporting conics are equal
    // up to a constant factor.
    Integer        factor1 = 1;
    Integer        factor2 = 1;

    if (CGAL::sign (this->_r) != ZERO)
      factor1 = this->_r;
    else if (CGAL::sign (this->_s) != ZERO)
      factor1 = this->_s;
    else if (CGAL::sign (this->_t) != ZERO)
      factor1 = this->_t;
    else if (CGAL::sign (this->_u) != ZERO)
      factor1 = this->_u;
    else if (CGAL::sign (this->_v) != ZERO)
      factor1 = this->_v;
    else if (CGAL::sign (this->_w) != ZERO)
      factor1 = this->_w;

    if (CGAL::sign (arc._r) != ZERO)
      factor2 = arc._r;
    else if (CGAL::sign (arc._s) != ZERO)
      factor2 = arc._s;
    else if (CGAL::sign (arc._t) != ZERO)
      factor2 = arc._t;
    else if (CGAL::sign (arc._u) != ZERO)

      factor2 = arc._u;
    else if (CGAL::sign (arc._v) != ZERO)
      factor2 = arc._v;
    else if (CGAL::sign (arc._w) != ZERO)
      factor2 = arc._w;

    return (CGAL::compare  (this->_r * factor2, arc._r * factor1) == EQUAL &&
            CGAL::compare  (this->_s * factor2, arc._s * factor1) == EQUAL &&
            CGAL::compare  (this->_t * factor2, arc._t * factor1) == EQUAL &&
            CGAL::compare  (this->_u * factor2, arc._u * factor1) == EQUAL &&
            CGAL::compare  (this->_v * factor2, arc._v * factor1) == EQUAL &&
            CGAL::compare  (this->_w * factor2, arc._w * factor1) == EQUAL);
  }

  /*!
   * Get the i'th order derivative by x of the conic at the point p=(x,y).
   * \param p The point where we derive.
   * \param i The order of the derivatives (either 1, 2 or 3).
   * \param slope_numer The numerator of the slope.
   * \param slope_denom The denominator of the slope.
   * \todo Allow higher order derivatives.
   */
  void _derive_by_x_at (const Point_2& p, const unsigned int& i,
                        Algebraic& slope_numer, Algebraic& slope_denom) const
  {
    if (_is_special_segment())
    {
      // Special treatment for special segments, given by (a*x + b*y + c = 0),
      // so their first-order derivative by x is simply -a/b. The higher-order
      // derivatives are all 0.
      if (i == 1)
      {
        if (CGAL::sign (this->_extra_data_P->b) != NEGATIVE)
        {
          slope_numer = - this->_extra_data_P->a;
          slope_denom = this->_extra_data_P->b;
        }
        else
        {
          slope_numer = this->_extra_data_P->a;
          slope_denom = - this->_extra_data_P->b;
        }
      }
      else
      {
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
    const Algebraic  _two = 2;
    const Algebraic  sl_numer = _two*alg_r*p.x() + alg_t*p.y() + alg_u;
    const Algebraic  sl_denom = _two*alg_s*p.y() + alg_t*p.x() + alg_v;

    if (i == 1)
    {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE)
      {
        slope_numer = -sl_numer;
        slope_denom = sl_denom;
      }
      else
      {
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
    const Algebraic  sl2_numer = alg_s * sl_numer*sl_numer -
                                 alg_t * sl_numer*sl_denom +
                                 alg_r * sl_denom*sl_denom;
    const Algebraic  sl2_denom = sl_denom*sl_denom*sl_denom;

    if (i == 2)
    {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE)
      {
        slope_numer = -_two *sl2_numer;
        slope_denom = sl2_denom;
      }
      else
      {
        slope_numer = _two *sl2_numer;
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
    const Algebraic  sl3_numer = (_two * alg_s * sl_numer -
                                  alg_t * sl_denom) * sl2_numer;
    const Algebraic  sl3_denom = sl_denom*sl_denom * sl2_denom;

    if (i == 3)
    {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE)
      {
        slope_numer = -6 * sl3_numer;
        slope_denom = sl3_denom;
      }
      else
      {
        slope_numer = 6 * sl3_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // \todo Handle higher-order derivatives as well.
    CGAL_error();
    return;
  }

  /*!
   * Get the i'th order derivative by y of the conic at the point p=(x,y).
   * \param p The point where we derive.
   * \param i The order of the derivatives (either 1, 2 or 3).
   * \param slope_numer The numerator of the slope.
   * \param slope_denom The denominator of the slope.
   * \todo Allow higher order derivatives.
   */
  void _derive_by_y_at (const Point_2& p, const int& i,
                        Algebraic& slope_numer, Algebraic& slope_denom) const
  {
    if (_is_special_segment())
    {
      // Special treatment for special segments, given by (a*x + b*y + c = 0),
      // so their first-order derivative by x is simply -b/a. The higher-order
      // derivatives are all 0.
      if (i == 1)
      {
        if (CGAL::sign (this->_extra_data_P->a) != NEGATIVE)
        {
          slope_numer = - this->_extra_data_P->b;
          slope_denom = this->_extra_data_P->a;
        }
        else
        {
          slope_numer = this->_extra_data_P->b;
          slope_denom = - this->_extra_data_P->a;
        }
      }
      else
      {
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
    const Algebraic  _two = 2;
    const Algebraic  sl_numer = _two*alg_s*p.y() + alg_t*p.x() + alg_v;
    const Algebraic  sl_denom = _two*alg_r*p.x() + alg_t*p.y() + alg_u;

    if (i == 1)
    {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE)
      {
        slope_numer = -sl_numer;
        slope_denom = sl_denom;
      }
      else
      {
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
    const Algebraic  sl2_numer = alg_r * sl_numer*sl_numer -
                                 alg_t * sl_numer*sl_denom +
                                 alg_s * sl_denom*sl_denom;
    const Algebraic  sl2_denom = sl_denom*sl_denom*sl_denom;

    if (i == 2)
    {

      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE)
      {
        slope_numer = -_two *sl2_numer;
        slope_denom = sl2_denom;
      }
      else
      {
        slope_numer = _two *sl2_numer;
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
    const Algebraic  sl3_numer = (_two * alg_r * sl_numer -
                                  alg_t * sl_denom) * sl2_numer;
    const Algebraic  sl3_denom = sl_denom*sl_denom * sl2_denom;

    if (i == 3)
    {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE)
      {
        slope_numer = -6 * sl3_numer;
        slope_denom = sl3_denom;
      }
      else
      {
        slope_numer = 6 * sl3_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // \todo Handle higher-order derivatives as well.
    CGAL_error();
    return;
  }

  /*!
   * Compute the overlap with a given arc, which is supposed to have the same
   * supporting conic curve as this arc.
   * \param arc The given arc.
   * \param overlap Output: The overlapping arc (if any).
   * \return Whether we found an overlap.
   */
  bool _compute_overlap (const Self& arc, Self& overlap) const
  {
    // Check if the two arcs are identical.
    if (equals (arc))
    {
      overlap = arc;
      return (true);
    }

    if (_is_strictly_between_endpoints (arc.left()))
    {
      if (_is_strictly_between_endpoints(arc.right()))
      {
        // Case 1 - *this:     +----------->
        //            arc:       +=====>
        overlap = arc;
        return (true);
      }
      else
      {
        // Case 2 - *this:     +----------->
        //            arc:               +=====>
        overlap = *this;

        if ((overlap._info & IS_DIRECTED_RIGHT) != 0)
          overlap._source = arc.left();
        else
          overlap._target = arc.left();

        return (true);
      }
    }
    else if (_is_strictly_between_endpoints (arc.right()))
    {
      // Case 3 - *this:     +----------->
      //            arc:   +=====>
      overlap = *this;

      if ((overlap._info & IS_DIRECTED_RIGHT) != 0)
        overlap._target = arc.right();
      else
        overlap._source = arc.right();

      return (true);
    }
    else if (arc._is_between_endpoints (this->_source) &&
             arc._is_between_endpoints (this->_target) &&
             (arc._is_strictly_between_endpoints(this->_source) ||
              arc._is_strictly_between_endpoints(this->_target)))
    {
      // Case 4 - *this:     +----------->
      //            arc:   +================>
      overlap = *this;
      return (true);
    }

    // If we reached here, there are no overlaps:
    return (false);
  }

  /*!
   * Intersect the supporing conic curves of this arc and the given arc.
   * \param arc The arc to intersect with.
   * \param inter_list The list of intersection points.
   */
  void _intersect_supporting_conics(const Self& arc,
                                    Intersection_list& inter_list) const
  {
    if (_is_special_segment() && ! arc._is_special_segment()) {
      // If one of the arcs is a special segment, make sure it is (arc).
      arc._intersect_supporting_conics(*this, inter_list);
      return;
    }

    const int deg1 = ((this->_info & DEGREE_MASK) == DEGREE_1) ? 1 : 2;
    const int deg2 = ((arc._info & DEGREE_MASK) == DEGREE_1) ? 1 : 2;
    Nt_traits nt_traits;
    Algebraic xs[4];
    int n_xs = 0;
    Algebraic ys[4];
    int n_ys = 0;

    if (arc._is_special_segment()) {
      // The second arc is a special segment (a*x + b*y + c = 0).
      if (_is_special_segment()) {
        // Both arc are sepcial segment, so they have at most one intersection
        // point.
        Algebraic denom = this->_extra_data_P->a * arc._extra_data_P->b -
          this->_extra_data_P->b * arc._extra_data_P->a;

        if (CGAL::sign (denom) != CGAL::ZERO) {
          xs[0] = (this->_extra_data_P->b * arc._extra_data_P->c -
                   this->_extra_data_P->c * arc._extra_data_P->b) / denom;
          n_xs = 1;

          ys[0] = (this->_extra_data_P->c * arc._extra_data_P->a -
                   this->_extra_data_P->a * arc._extra_data_P->c) / denom;
          n_ys = 1;
        }
      }
      else {
        // Compute the x-coordinates of the intersection points.
        n_xs = _compute_resultant_roots (nt_traits,
                                         alg_r, alg_s, alg_t,
                                         alg_u, alg_v, alg_w,
                                         deg1,
                                         arc._extra_data_P->a,
                                         arc._extra_data_P->b,
                                         arc._extra_data_P->c,
                                         xs);
        CGAL_assertion (n_xs <= 2);

        // Compute the y-coordinates of the intersection points.
        n_ys = _compute_resultant_roots (nt_traits,
                                         alg_s, alg_r, alg_t,
                                         alg_v, alg_u, alg_w,
                                         deg1,
                                         arc._extra_data_P->b,
                                         arc._extra_data_P->a,
                                         arc._extra_data_P->c,
                                         ys);
        CGAL_assertion (n_ys <= 2);
      }
    }
    else {
      // Compute the x-coordinates of the intersection points.
      n_xs = _compute_resultant_roots (nt_traits,
                                       this->_r, this->_s, this->_t,
                                       this->_u, this->_v, this->_w,
                                       deg1,
                                       arc._r, arc._s, arc._t,
                                       arc._u, arc._v, arc._w,
                                       deg2,
                                       xs);
      CGAL_assertion (n_xs <= 4);

      // Compute the y-coordinates of the intersection points.
      n_ys = _compute_resultant_roots (nt_traits,
                                       this->_s, this->_r, this->_t,
                                       this->_v, this->_u, this->_w,
                                       deg1,
                                       arc._s, arc._r, arc._t,
                                       arc._v, arc._u, arc._w,
                                       deg2,
                                       ys);
      CGAL_assertion(n_ys <= 4);
    }

    // Pair the coordinates of the intersection points. As the vectors of
    // x and y-coordinates are sorted in ascending order, we output the
    // intersection points in lexicographically ascending order.
    unsigned int  mult;
    int   i, j;

    if (arc._is_special_segment())
    {
      if (n_xs == 0 || n_ys == 0)
        return;

      if (n_xs == 1 && n_ys == 1)
      {
        // Single intersection.
        Conic_point_2 ip (xs[0], ys[0]);

        ip.set_generating_conic (_id);
        ip.set_generating_conic (arc._id);

        // In case the other curve is of degree 2, this is a tangency point.
        mult = (deg1 == 1 || _is_special_segment()) ? 1 : 2;
        inter_list.push_back(Intersection_point (ip, mult));
      }
      else if (n_xs == 1 && n_ys == 2)
      {
        Conic_point_2         ip1 (xs[0], ys[0]);

        ip1.set_generating_conic (_id);
        ip1.set_generating_conic (arc._id);

        inter_list.push_back(Intersection_point (ip1, 1));

        Conic_point_2         ip2 (xs[0], ys[1]);

        ip2.set_generating_conic (_id);
        ip2.set_generating_conic (arc._id);

        inter_list.push_back(Intersection_point (ip2, 1));
      }
      else if (n_xs == 2 && n_ys == 1)
      {
        Conic_point_2         ip1 (xs[0], ys[0]);

        ip1.set_generating_conic (_id);
        ip1.set_generating_conic (arc._id);

        inter_list.push_back(Intersection_point (ip1, 1));

        Conic_point_2         ip2 (xs[1], ys[0]);

        ip2.set_generating_conic (_id);
        ip2.set_generating_conic (arc._id);

        inter_list.push_back(Intersection_point (ip2, 1));

      }
      else {
        CGAL_assertion (n_xs == 2 && n_ys == 2);

        // The x-coordinates and the y-coordinates are given in ascending
        // order. If the slope of the segment is positive, we pair the
        // coordinates as is - otherwise, we swap the pairs.
        int ind_first_y = 0, ind_second_y = 1;

        if (CGAL::sign (arc._extra_data_P->b) ==
            CGAL::sign(arc._extra_data_P->a))
        {
          ind_first_y = 1;
          ind_second_y = 0;
        }

        Conic_point_2 ip1(xs[0], ys[ind_first_y]);

        ip1.set_generating_conic(_id);
        ip1.set_generating_conic(arc._id);

        inter_list.push_back(Intersection_point (ip1, 1));

        Conic_point_2 ip2(xs[1], ys[ind_second_y]);

        ip2.set_generating_conic (_id);
        ip2.set_generating_conic (arc._id);

        inter_list.push_back(Intersection_point(ip2, 1));
      }

      return;
    }

    for (i = 0; i < n_xs; i++) {
      for (j = 0; j < n_ys; j++) {
        if (_is_on_supporting_conic (xs[i], ys[j]) &&
            arc._is_on_supporting_conic (xs[i], ys[j]))
        {
          // Create the intersection point and set its generating conics.
          Conic_point_2 ip(xs[i], ys[j]);

          ip.set_generating_conic (_id);
          ip.set_generating_conic (arc._id);

          // Compute the multiplicity of the intersection point.
          if (deg1 == 1 && deg2 == 1) mult = 1;
          else mult = _multiplicity_of_intersection_point(arc, ip);

          // Insert the intersection point to the output list.
          inter_list.push_back(Intersection_point(ip, mult));
        }
      }
    }

    return;
  }

  /*!
   * Compute the multiplicity of an intersection point.
   * \param arc The arc to intersect with.
   * \param p The intersection point.
   * \return The multiplicity of the intersection point.
   */
  unsigned int _multiplicity_of_intersection_point (const Self& arc,
                                                    const Point_2& p) const
  {
    CGAL_assertion (! _is_special_segment() || ! arc._is_special_segment());

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic      slope1_numer, slope1_denom;
    Algebraic      slope2_numer, slope2_denom;

    _derive_by_x_at (p, 1, slope1_numer, slope1_denom);
    arc._derive_by_x_at (p, 1, slope2_numer, slope2_denom);

    if (CGAL::compare (slope1_numer*slope2_denom,
                       slope2_numer*slope1_denom) != EQUAL)
    {
      // Different slopes at p - the mutiplicity of p is 1:
      return (1);
    }

    if (CGAL::sign (slope1_denom) != ZERO &&
        CGAL::sign (slope2_denom) != ZERO)
    {
      // The curves do not have a vertical slope at p.
      // Compare their second-order derivative by x:
      _derive_by_x_at (p, 2, slope1_numer, slope1_denom);
      arc._derive_by_x_at (p, 2, slope2_numer, slope2_denom);
    }
    else
    {
      // Both curves have a vertical slope at p.
      // Compare their second-order derivative by y:
      _derive_by_y_at (p, 2, slope1_numer, slope1_denom);
      arc._derive_by_y_at (p, 2, slope2_numer, slope2_denom);
    }

    if (CGAL::compare (slope1_numer*slope2_denom,
                       slope2_numer*slope1_denom) != EQUAL)
    {
      // Different curvatures at p - the mutiplicity of p is 2:
      return (2);
    }

    // If we reached here, the multiplicity of the intersection point is 3:
    return (3);
  }
  //@}

};

/*!
 * Exporter for x-monotone conic arcs.
 */
template <class Conic_arc_2>
std::ostream& operator<< (std::ostream& os,
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

  if (arc.orientation() == CLOCKWISE)
    os << "--cw-->";
  else if (arc.orientation() == COUNTERCLOCKWISE)
    os << "--ccw-->";
  else
    os << "--l-->";

  os << " (" << CGAL::to_double(arc.target().x()) << ","
     << CGAL::to_double(arc.target().y()) << ")";

  return (os);
}

} //namespace CGAL

#endif
