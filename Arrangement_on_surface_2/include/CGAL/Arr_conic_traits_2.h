// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein   <wein@post.tau.ac.il>
//            Waqar Khan <wkhan@mpi-inf.mpg.de>
//            Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_ARR_CONIC_TRAITS_2_H
#define CGAL_ARR_CONIC_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The conic traits-class for the arrangement package.
 */

// Keep the following 2 lines first.
#include <cmath>
#include <fstream>
#include <atomic>
#include <memory>
#include <map>

#include <boost/math/constants/constants.hpp>

#include <CGAL/Cartesian.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_geometry_traits/Conic_arc_2.h>
#include <CGAL/Arr_geometry_traits/Conic_x_monotone_arc_2.h>
#include <CGAL/Arr_geometry_traits/Conic_point_2.h>
#include <CGAL/Arr_geometry_traits/Conic_intersections_2.h>
#include <CGAL/Bbox_2.h>

namespace CGAL {

/*! \class A traits class for maintaining an arrangement of conic arcs (bounded
 * segments of algebraic curves of degree 2 at most).
 *
 * The class is templated with two parameters:
 * Rat_kernel A kernel that provides the input objects or coefficients.
 *            Rat_kernel::FT should be an integral or a rational type.
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices, which are algebraic
 *            numbers of degree up to 4 (preferably it is CORE::Expr).
 * Nt_traits A traits class for performing various operations on the integer,
 *           rational and algebraic types.
 */
template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Arr_conic_traits_2 {
public:
  typedef RatKernel                       Rat_kernel;
  typedef AlgKernel                       Alg_kernel;
  typedef NtTraits                        Nt_traits;

  typedef typename Rat_kernel::FT         Rational;
  typedef typename Rat_kernel::Point_2    Rat_point_2;
  typedef typename Rat_kernel::Segment_2  Rat_segment_2;
  typedef typename Rat_kernel::Line_2     Rat_line_2;
  typedef typename Rat_kernel::Circle_2   Rat_circle_2;

  typedef typename Alg_kernel::FT         Algebraic;
  typedef typename Alg_kernel::Point_2    Alg_point_2;

  typedef typename Nt_traits::Integer     Integer;

  // Category tags:
  typedef Tag_true                        Has_left_category;
  typedef Tag_true                        Has_merge_category;
  typedef Tag_false                       Has_do_intersect_category;
  //typedef std::true_type                Has_line_segment_constructor;

  typedef Arr_oblivious_side_tag          Left_side_category;
  typedef Arr_oblivious_side_tag          Bottom_side_category;
  typedef Arr_oblivious_side_tag          Top_side_category;
  typedef Arr_oblivious_side_tag          Right_side_category;

  // Traits objects:
  typedef Conic_arc_2<Rat_kernel, Alg_kernel, Nt_traits>  Curve_2;
  typedef Conic_x_monotone_arc_2<Curve_2>                 X_monotone_curve_2;
  typedef Conic_point_2<Alg_kernel>                       Point_2;
  typedef size_t                                          Multiplicity;

private:
  // Type definition for the intersection points mapping.
  using Conic_id = typename Point_2::Conic_id;
  using Conic_pair = std::pair<Conic_id, Conic_id>;

  /*! \struct Less functor for Conic_pair.
   */
  struct Less_conic_pair {
    bool operator()(const Conic_pair& cp1, const Conic_pair& cp2) const {
      // Compare the pairs of IDs lexicographically.
      return ((cp1.first < cp2.first) ||
              ((cp1.first == cp2.first) && (cp1.second < cp2.second)));
    }
  };

  typedef std::pair<Point_2, Multiplicity>          Intersection_point;
  typedef std::list<Intersection_point>             Intersection_list;
  typedef std::map<Conic_pair, Intersection_list, Less_conic_pair>
                                                    Intersection_map;
  typedef typename Intersection_map::iterator       Intersection_map_iterator;


  typedef std::shared_ptr<Rat_kernel>               Shared_rat_kernel;
  typedef std::shared_ptr<Alg_kernel>               Shared_alg_kernel;
  typedef std::shared_ptr<Nt_traits>                Shared_nt_traits;

  const Shared_rat_kernel m_rat_kernel;
  const Shared_alg_kernel m_alg_kernel;
  const Shared_nt_traits m_nt_traits;

  mutable Intersection_map m_inter_map; // Mapping conic pairs to their
                                        // intersection points.

public:
  /*! Default constructor.
   */
  Arr_conic_traits_2() {}

  /*! Construct from resources.
   */
  Arr_conic_traits_2(Shared_rat_kernel rat_kernel,
                     Shared_alg_kernel alg_kernel,
                     Shared_nt_traits nt_traits) :
    m_rat_kernel(rat_kernel),
    m_alg_kernel(alg_kernel),
    m_nt_traits(nt_traits)
  {}

  /*! Obtain the rational kernel.
   */
  Shared_rat_kernel rat_kernel() const { return m_rat_kernel; }

  /*! Obtain the algebraic kernel.
   */
  Shared_alg_kernel alg_kernel() const { return m_alg_kernel; }

  /*! Obtain the nt traits.
   */
  Shared_nt_traits nt_traits() const { return m_nt_traits; }

  /*! Obtain the next conic index. */
  static size_t get_index() {
#ifdef CGAL_NO_ATOMIC
    static size_t index;
#else
    static std::atomic<size_t> index;
#endif
    return (++index);
  }

  /// \name Basic functor definitions.
  //@{

  class Compare_x_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits.
     */
    Compare_x_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) < x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    { return m_traits.m_alg_kernel->compare_x_2_object()(p1, p2); }
  };

  /*! Obtain a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(*this); }

  class Compare_xy_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits.
     */
    Compare_xy_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Compares two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return `LARGER` if `x(p1) > x(p2)`, or if `x(p1) = x(p2)` and `y(p1) > y(p2)`;
     *         `SMALLER` if `x(p1) < x(p2)`, or if `x(p1) = x(p2)` and `y(p1) < y(p2)`;
     *         `EQUAL` if the two points are equal.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { return m_traits.m_alg_kernel->compare_xy_2_object()(p1, p2); }
  };

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(*this); }

  class Construct_min_vertex_2 {
  public:
    /*! Obtain the left endpoint of an x-monotone arc.
     * \param cv The arc.
     * \return The left endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& xcv) const
    { return xcv.left(); }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  class Construct_max_vertex_2 {
  public:
    /*! Obtain the right endpoint of the x-monotone arc.
     * \param cv The arc.
     * \return The right endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& xcv) const
    { return xcv.right(); }
  };

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  class Is_vertical_2 {
  public:
    /*! Check whether a given x-monotone arc is a vertical segment.
     * \param cv The vertical segment.
     * \return `true` if the arc is a vertical segment; `false` otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv) const
    { return cv.is_vertical(); }
  };

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(); }

  class Compare_y_at_x_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits.
     */
    Compare_y_at_x_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Return the location of a given point with respect to an input arc.
     * \param xcv The arc.
     * \param p The point.
     * \pre `p` is in the \f$x\f$-range of `xcv`.
     * \return `SMALLER` if `y(p) < xcv(x(p))`, i.e. the point is below the arc;
     *         `LARGER` if `y(p) > xcv(x(p))`, i.e. the point is above the arc;
     *         `EQUAL` if `p` lies on the curve.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv) const {
      auto cmp_y = m_traits.m_alg_kernel->compare_y_2_object();

      if (xcv.is_vertical()) {
        // A special treatment for vertical segments:
        // In case p has the same x c-ordinate of the vertical segment, compare
        // it to the segment endpoints to determine its position.
        Comparison_result res1 = cmp_y(p, xcv.left());
        Comparison_result res2 = cmp_y(p, xcv.right());
        return (res1 == res2) ? res1 : EQUAL;
      }

      // Check whether the point is exactly on the curve.
      if (m_traits.contains_point(xcv, p)) return EQUAL;

      // Obtain a point q on the x-monotone arc with the same x coordinate as p.
      Point_2 q;

      auto cmp_x = m_traits.m_alg_kernel->compare_x_2_object();
      Comparison_result x_res_left = cmp_x(p, xcv.left());
      if (x_res_left == EQUAL) q = xcv.left();
      else {
        CGAL_precondition(x_res_left != SMALLER);
        auto x_res_right = cmp_x(p, xcv.right());
        if (x_res_right == EQUAL) q = xcv.right();
        else {
          CGAL_precondition(x_res_right != LARGER);
          q = m_traits.point_at_x(xcv, p);
        }
      }

      // Compare p with the a point of the curve with the same x coordinate.
      return cmp_y(p, q);
    }
  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(*this); }

  class Compare_y_at_x_left_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits.
     */
    Compare_y_at_x_left_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Compares the \f$y\f$ value of two \f$x\f$-monotone arcs immediately
     * to the left of their intersection point.
     * \param xcv1 The first arc.
     * \param xcv2 The second arc.
     * \param p The intersection point.
     * \pre The point `p` lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of `xcv1` with respect to `xcv2` immdiately
     *         to the left of `p`: `SMALLER`, `LARGER`, or `EQUAL`.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 const Point_2& p) const
    {
      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition(m_traits.contains_point(xcv1, p) &&
                        m_traits.contains_point(xcv2, p));

      CGAL_precondition_code(const auto ker = m_traits.m_alg_kernel);
      CGAL_precondition(ker->compare_xy_2_object()(p, xcv1.left()) == LARGER &&
                        ker->compare_xy_2_object()(p, xcv2.left()) == LARGER);

      // If one of the curves is vertical, it is below the other one.
      if (xcv1.is_vertical()) return (xcv2.is_vertical()) ? EQUAL : SMALLER;
      else if (xcv2.is_vertical()) return LARGER;

      // Compare the two curves immediately to the left of p:
      return compare_to_left(xcv1, xcv2, p);
    }

  private:
    /*! Compare two arcs immediately to the leftt of their intersection point.
     * \param xcv1 The first compared arc.
     * \param xcv2 The second compared arc.
     * \param p The reference intersection point.
     * \return The relative position of the arcs to the left of `p`.
     * \pre Both arcs we compare are not vertical segments.
     */
    Comparison_result compare_to_left(const X_monotone_curve_2& xcv1,
                                      const X_monotone_curve_2& xcv2,
                                      const Point_2& p) const {
      CGAL_precondition(! xcv1.is_vertical() && ! xcv2.is_vertical());

      // In case one arc is facing upwards and another facing downwards, it is
      // clear that the one facing upward is above the one facing downwards.
      if (m_traits.has_same_supporting_conic(xcv1, xcv2)) {
        if (xcv1.test_flag(X_monotone_curve_2::FACING_UP) &&
            xcv2.test_flag(X_monotone_curve_2::FACING_DOWN))
          return LARGER;
        else if (xcv1.test_flag(X_monotone_curve_2::FACING_DOWN) &&
                 xcv2.test_flag(X_monotone_curve_2::FACING_UP))
          return SMALLER;

        // In this case the two arcs overlap.
        CGAL_assertion(xcv1.facing_mask() == xcv2.facing_mask());

        return EQUAL;
      }

      // Compare the slopes of the two arcs at p, using their first-order
      // partial derivatives.
      Algebraic slope1_numer, slope1_denom;
      Algebraic slope2_numer, slope2_denom;

      xcv1.derive_by_x_at(p, 1, slope1_numer, slope1_denom);
      xcv2.derive_by_x_at(p, 1, slope2_numer, slope2_denom);

      // Check if any of the slopes is vertical.
      const bool is_vertical_slope1 = (CGAL::sign (slope1_denom) == ZERO);
      const bool is_vertical_slope2 = (CGAL::sign (slope2_denom) == ZERO);

      if (! is_vertical_slope1 && ! is_vertical_slope2) {
        // The two derivatives at p are well-defined: use them to determine
        // which arc is above the other (the one with a larger slope is below).
        Comparison_result slope_res = CGAL::compare(slope2_numer*slope1_denom,
                                                    slope1_numer*slope2_denom);

        if (slope_res != EQUAL) return slope_res;

        // Use the second-order derivative.
        xcv1.derive_by_x_at(p, 2, slope1_numer, slope1_denom);
        xcv2.derive_by_x_at(p, 2, slope2_numer, slope2_denom);

        slope_res = CGAL::compare(slope1_numer*slope2_denom,
                                  slope2_numer*slope1_denom);

        if (slope_res != EQUAL) return (slope_res);

        // Use the third-order derivative.
        xcv1.derive_by_x_at(p, 3, slope1_numer, slope1_denom);
        xcv2.derive_by_x_at(p, 3, slope2_numer, slope2_denom);

        slope_res = CGAL::compare(slope2_numer*slope1_denom,
                                  slope1_numer*slope2_denom);

        // \todo Handle higher-order derivatives:
        CGAL_assertion(slope_res != EQUAL);

        return slope_res;
      }
      else if (! is_vertical_slope2) {
        // The first arc has a vertical slope at p: check whether it is
        // facing upwards or downwards and decide accordingly.
        CGAL_assertion(xcv1.facing_mask() != 0);

        return (xcv1.test_flag(X_monotone_curve_2::FACING_UP)) ?
          LARGER : SMALLER;
      }
      else if (! is_vertical_slope1) {
        // The second arc has a vertical slope at p_int: check whether it is
        // facing upwards or downwards and decide accordingly.
        CGAL_assertion(xcv2.facing_mask() != 0);

        return (xcv2.test_flag(X_monotone_curve_2::FACING_UP)) ?
          SMALLER : LARGER;
      }

      // The two arcs have vertical slopes at p_int:
      // First check whether one is facing up and one down. In this case the
      // comparison result is trivial.
      if (xcv1.test_flag(X_monotone_curve_2::FACING_UP) &&
          xcv2.test_flag(X_monotone_curve_2::FACING_DOWN))
        return LARGER;
      else if (xcv1.test_flag(X_monotone_curve_2::FACING_DOWN) &&
               xcv2.test_flag(X_monotone_curve_2::FACING_UP))
        return SMALLER;

      // Compute the second-order derivative by y and act according to it.
      xcv1.derive_by_y_at(p, 2, slope1_numer, slope1_denom);
      xcv2.derive_by_y_at(p, 2, slope2_numer, slope2_denom);

      Comparison_result slope_res =
        CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);

      // If necessary, use the third-order derivative by y.
      if (slope_res == EQUAL) {
        // \todo Check this!
        xcv1.derive_by_y_at(p, 3, slope1_numer, slope1_denom);
        xcv2.derive_by_y_at(p, 3, slope2_numer, slope2_denom);

        slope_res =
          CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);
      }

      // \todo Handle higher-order derivatives:
      CGAL_assertion(slope_res != EQUAL);

      // Check whether both are facing up.
      if (xcv1.test_flag(X_monotone_curve_2::FACING_UP) &&
          xcv2.test_flag(X_monotone_curve_2::FACING_UP))
        return ((slope_res == LARGER) ? SMALLER : LARGER);

      // Both are facing down.
      return slope_res;
    }

  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(*this); }

  class Compare_y_at_x_right_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits.
     */
    Compare_y_at_x_right_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Compares the `y` value of two \f$x\f$-monotone arcs immediately
     * to the right of their intersection point.
     * \param xcv1 The first arc.
     * \param xcv2 The second arc.
     * \param p The intersection point.
     * \pre The point `p` lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of `xcv1` with respect to `xcv2` immdiately
     *         to the right of `p`: `SMALLER`, `LARGER`, or `EQUAL`.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 const Point_2& p) const
    {
      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition(m_traits.contains_point(xcv1, p) &&
                        m_traits.contains_point(xcv2, p));
      CGAL_precondition_code(const auto ker = m_traits.m_alg_kernel);
      CGAL_precondition_code(auto cmp_xy = ker->compare_xy_2_object());
      CGAL_precondition(cmp_xy(p, xcv1.right()) == SMALLER &&
                        cmp_xy(p, xcv2.right()) == SMALLER);

      // If one of the curves is vertical, it is above the other one.
      if (xcv1.is_vertical()) return (xcv2.is_vertical()) ? EQUAL : LARGER;
      else if (xcv2.is_vertical()) return SMALLER;

      // Compare the two curves immediately to the right of p:
      return compare_to_right(xcv1, xcv2, p);
    }

  private:
  /*! Compare two arcs immediately to the right of their intersection point.
   * \param xcv1 The first compared arc.
   * \param xcv2 The second compared arc.
   * \param p The reference intersection point.
   * \return The relative position of the arcs to the right of `p`.
   * \pre Both arcs we compare are not vertical segments.
   */
  Comparison_result compare_to_right(const X_monotone_curve_2& xcv1,
                                     const X_monotone_curve_2& xcv2,
                                     const Point_2& p) const {
    CGAL_precondition(! xcv1.is_vertical() && ! xcv2.is_vertical());

    // In case one arc is facing upwards and another facing downwards, it is
    // clear that the one facing upward is above the one facing downwards.
    if (m_traits.has_same_supporting_conic(xcv1, xcv2)) {
      if (xcv1.test_flag(X_monotone_curve_2::FACING_UP) &&
          xcv2.test_flag(X_monotone_curve_2::FACING_DOWN))
        return LARGER;
      else if (xcv1.test_flag(X_monotone_curve_2::FACING_DOWN) &&
               xcv2.test_flag(X_monotone_curve_2::FACING_UP))
        return SMALLER;

      // In this case the two arcs overlap.
      CGAL_assertion(xcv1.facing_mask() == xcv2.facing_mask());
      return EQUAL;
    }

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic slope1_numer, slope1_denom;
    Algebraic slope2_numer, slope2_denom;

    xcv1.derive_by_x_at(p, 1, slope1_numer, slope1_denom);
    xcv2.derive_by_x_at(p, 1, slope2_numer, slope2_denom);

    // Check if any of the slopes is vertical.
    const bool is_vertical_slope1 = (CGAL::sign(slope1_denom) == ZERO);
    const bool is_vertical_slope2 = (CGAL::sign(slope2_denom) == ZERO);

    if (! is_vertical_slope1 && ! is_vertical_slope2) {
      // The two derivatives at p are well-defined: use them to determine
      // which arc is above the other (the one with a larger slope is below).
      Comparison_result slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the second-order derivative.
      xcv1.derive_by_x_at(p, 2, slope1_numer, slope1_denom);
      xcv2.derive_by_x_at(p, 2, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the third-order derivative.
      xcv1.derive_by_x_at(p, 3, slope1_numer, slope1_denom);
      xcv2.derive_by_x_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      // \todo Handle higher-order derivatives:
      CGAL_assertion(slope_res != EQUAL);

      return slope_res;
    }
    else if (! is_vertical_slope2) {
      // The first arc has a vertical slope at p: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion(xcv1.facing_mask() != 0);

      return (xcv1.test_flag(X_monotone_curve_2::FACING_UP)) ? LARGER : SMALLER;
    }
    else if (! is_vertical_slope1) {
      // The second arc has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion(xcv2.facing_mask() != 0);

      return (xcv2.test_flag(X_monotone_curve_2::FACING_UP)) ? SMALLER : LARGER;
    }

    // The two arcs have vertical slopes at p_int:
    // First check whether one is facing up and one down. In this case the
    // comparison result is trivial.
    if (xcv1.test_flag(X_monotone_curve_2::FACING_UP) &&
        xcv2.test_flag(X_monotone_curve_2::FACING_DOWN)) return LARGER;
    else if (xcv1.test_flag(X_monotone_curve_2::FACING_DOWN) &&
             xcv2.test_flag(X_monotone_curve_2::FACING_UP)) return SMALLER;

    // Compute the second-order derivative by y and act according to it.
    xcv1.derive_by_y_at(p, 2, slope1_numer, slope1_denom);
    xcv2.derive_by_y_at(p, 2, slope2_numer, slope2_denom);

    Comparison_result slope_res =
      CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

    // If necessary, use the third-order derivative by y.
    if (slope_res == EQUAL) {
      // \todo Check this!
      xcv1.derive_by_y_at(p, 3, slope1_numer, slope1_denom);
      xcv2.derive_by_y_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);
    }

    // \todo Handle higher-order derivatives:
    CGAL_assertion(slope_res != EQUAL);

    if (xcv1.test_flag(X_monotone_curve_2::FACING_UP) &&
        xcv2.test_flag(X_monotone_curve_2::FACING_UP))
      return (slope_res == LARGER) ? SMALLER : LARGER;  // both are facing up
    return slope_res;                                   // both are facing down
  }

  };

  /*! Obtain a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(*this); }

  class Equal_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits.
     */
    Equal_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Check whether two \f$x\f$-monotone curves are the same (have the same
     * graph).
     * \param xcv1 The first curve.
     * \param xcv2 The second curve.
     * \return `true` if the two curves are the same; `false` otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    {
      if (&xcv1 == &xcv2) return true;
      return equals(xcv1, xcv2);
    }

    /*! Check whether two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return `true` if the two point are the same; `false` otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const {
      if (&p1 == &p2) return (true);
      return(m_traits.m_alg_kernel->compare_xy_2_object()(p1, p2) == EQUAL);
    }

  private:
    /*! Check whether two arcs are equal (have the same graph).
     * \param xcv1 The first compared arc.
     * \param xcv2 The second compared arc.
     * \return `true` if the two arcs have the same graph; `false` otherwise.
     */
    bool equals(const X_monotone_curve_2& xcv1,
                const X_monotone_curve_2& xcv2) const {
      // The two arc must have the same supporting conic curves.
      if (! m_traits.has_same_supporting_conic(xcv1, xcv2)) return false;

      auto eq = m_traits.m_alg_kernel->equal_2_object();

      // Check that the arc endpoints are the same.
      if (xcv1.orientation() == COLLINEAR) {
        CGAL_assertion(xcv2.orientation() == COLLINEAR);
        return((eq(xcv1.source(), xcv2.source()) &&
                eq(xcv1.target(), xcv2.target())) ||
               (eq(xcv1.source(), xcv2.target()) &&
                eq(xcv1.target(), xcv2.source())));
      }

      if (xcv1.orientation() == xcv2.m_orient) {
        // Same orientation - the source and target points must be the same.
        return (eq(xcv1.source(), xcv2.source()) &&
                eq(xcv1.target(), xcv2.target()));
      }

      // Reverse orientation - the source and target points must be swapped.
      return (eq(xcv1.source(), xcv2.target()) &&
              eq(xcv1.target(), xcv2.source()));
    }
  };

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object() const { return Equal_2(*this); }
  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! A function object that obtains the parameter space of a geometric
   * entity along the \f$x-\f$xaxis.
   */
  class Parameter_space_in_x_2 {
  public:
    /*! Obtains the parameter space at the end of an arc along the \f$x\f$-axis.
     * \param xcv The arc.
     * \param ce The arc end indicator:
     *    `ARR_MIN_END`&mdash;the minimal end of `xcv` or
     *    `ARR_MAX_END`&mdash;the maximal end of `xcv`.
     * \return the parameter space at the `ce` end of the arc `xcv`.
     *   `ARR_LEFT_BOUNDARY` &mdash;the arc approaches the identification curve from
     *                              the right at the arc left end.
     *   `ARR_INTERIOR`      &mdash;the arc does not approache the identification curve.
     *   `ARR_RIGHT_BOUNDARY`&mdash;the arc approaches the identification curve from
     *                              the left at the arc right end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & /* xcv */,
                                   Arr_curve_end /* ce */) const {
      CGAL_error_msg("Not implemented yet!");
      return ARR_INTERIOR;
    }

    /*! Obtains the parameter space at a point along the \f$x\f$-axis.
     * \param p The point.
     * \return the parameter space at `p`.
     */
    Arr_parameter_space operator()(const Point_2 ) const
    { return ARR_INTERIOR; }
  };

  /*! Obtain a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(); }

  /*! A function object that obtains the parameter space of a geometric
   * entity along the y-axis
   */
  class Parameter_space_in_y_2 {
  public:
    /*! Obtains the parameter space at the end of an arc along the \f$y\f$-axis .
     * Note that if the arc end coincides with a pole, then unless the arc
     * coincides with the identification curve, the arc end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the arc coincides with the identification curve, it is assumed to
     * be smaller than any other object.
     * \param xcv The arc.
     * \param ce The arc end indicator:
     *   `ARR_MIN_END`&mdash;the minimal end of `xcv` or
     *   `ARR_MAX_END`&mdash;the maximal end of `xcv`.
     * \return the parameter space at the `ce` end of the arc `xcv`.
     *   `ARR_BOTTOM_BOUNDARY`&mdash;the arc approaches the south pole at the arc
     *                               left end.
     *   `ARR_INTERIOR`       &mdash;the arc does not approache a contraction point.
     *   `ARR_TOP_BOUNDARY`   &mdash;the arc approaches the north pole at the arc
     *                               right end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& /* xcv */,
                                   Arr_curve_end /* ce */) const {
      CGAL_error_msg("Not implemented yet!");
      return ARR_INTERIOR;
    }

    /*! Obtains the parameter space at a point along the \f$y\f$-axis.
     * \param p The point.
     * \return The parameter space at `p`.
     */
    Arr_parameter_space operator()(const Point_2 /* p */) const
    { return ARR_INTERIOR; }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }

  //@}

  /// \name Intersections, subdivisions, and mergings
  //@{

  /*! \class Make_x_monotone_2
   * A functor for subdividing curves into x-monotone curves.
   */
  class Make_x_monotone_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits.
     */
    Make_x_monotone_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Subdivide a given conic arc into \f$x\f$-monotone sub arcs
     * and insert them to a given output iterator.
     * \param cv The arc.
     * \param oi the output iterator for the result. Its dereference type is a
     *           variant that wraps a \c Point_2 or an \c X_monotone_curve_2
     *           objects.
     * \return the past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {

      auto ctr_xcv = m_traits.construct_x_monotone_curve_2_object();

      // Increment the serial number of the curve cv, which will serve as its
      // unique identifier.
      auto index = Traits::get_index();
      Conic_id conic_id(index);

      // Find the points of vertical tangency to cv and act accordingly.
      Alg_point_2 vtan_ps[2];
      auto n_vtan_ps = m_traits.vertical_tangency_points(cv, vtan_ps);
      if (n_vtan_ps == 0) {
        // In case the given curve is already x-monotone:
        *oi++ = ctr_xcv(cv, conic_id);
        return oi;
      }

      // Split the conic arc into x-monotone sub-curves.
      if (cv.is_full_conic()) {
        // Make sure we have two vertical tangency points.
        CGAL_assertion(n_vtan_ps == 2);

        // In case the curve is a full conic, split it into two x-monotone
        // arcs, one going from ps[0] to ps[1], and the other from ps[1] to
        // ps[0].
        *oi++ = ctr_xcv(cv, vtan_ps[0], vtan_ps[1], conic_id);
        *oi++ = ctr_xcv(cv, vtan_ps[1], vtan_ps[0], conic_id);
      }
      else {
        if (n_vtan_ps == 1) {
          // Split the arc into two x-monotone sub-curves: one going from the
          // arc source to ps[0], and the other from ps[0] to the target.
          *oi++ = ctr_xcv(cv, cv.source(), vtan_ps[0], conic_id);
          *oi++ = ctr_xcv(cv, vtan_ps[0], cv.target(), conic_id);
        }
        else {
          CGAL_assertion(n_vtan_ps == 2);

          // Identify the first point we encounter when going from cv's source
          // to its target, and the second point we encounter. Note that the
          // two endpoints must both be below the line connecting the two
          // tangnecy points (or both lies above it).
          int ind_first = 0;
          int ind_second = 1;
          auto ker = m_traits.m_alg_kernel;
          auto cmp_y_at_x_2 = ker->compare_y_at_x_2_object();
          auto line = ker->construct_line_2_object()(vtan_ps[0], vtan_ps[1]);
          auto start_pos = cmp_y_at_x_2(cv.source(), line);
          auto order_vpts = ker->compare_x_2_object()(vtan_ps[0], vtan_ps[1]);

          CGAL_assertion((start_pos != EQUAL) &&
                         (cmp_y_at_x_2(cv.target(), line) == start_pos));
          CGAL_assertion(order_vpts != EQUAL);

          if (((cv.orientation() == COUNTERCLOCKWISE) &&
               (start_pos == order_vpts)) ||
              ((cv.orientation() == CLOCKWISE) && (start_pos != order_vpts)))
          {
            ind_first = 1;
            ind_second = 0;
          }

          // Split the arc into three x-monotone sub-curves.
          *oi++ = ctr_xcv(cv, cv.source(),
                          vtan_ps[ind_first],
                          conic_id);

          *oi++ = ctr_xcv(cv, vtan_ps[ind_first],
                          vtan_ps[ind_second],
                          conic_id);

          *oi++ = ctr_xcv(cv, vtan_ps[ind_second],
                          cv.target(), conic_id);
        }
      }

      return oi;
    }
  };

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(*this); }

  class Split_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits.
     */
    Split_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Split a given \f$x\f$-monotone arc at a given point into two sub-arcs.
     * \param xcv The arc to split
     * \param p The split point.
     * \param xcv1 Output: The left resulting sub-arc (`p` is its right endpoint).
     * \param xcv2 Output: The right resulting sub-arc (`p` is its left endpoint).
     * \pre `p` lies on `xcv` but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& xcv, const Point_2 & p,
                    X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const
    { split(xcv, p, xcv1, xcv2); }

  private:
    /*! Split the arc into two at a given split point.
     * \param p The split point.
     * \param xcv1 Output: The first resulting arc, lying to the left of `p`.
     * \param xcv2 Output: The first resulting arc, lying to the right of `p`.
     * \pre `p` lies in the interior of the arc (not one of its endpoints).
     */
    void split(const X_monotone_curve_2& xcv, const Point_2& p,
               X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const {
      // Make sure that p lies on the interior of the arc.
      CGAL_precondition_code(auto eq = m_traits.m_alg_kernel->equal_2_object());
      CGAL_precondition(m_traits.contains_point(xcv, p) &&
                        ! eq(p, xcv.source()) && ! eq(p, xcv.target()));

      // Make copies of the current arc.
      xcv1 = xcv;
      xcv2 = xcv;

      // Assign the endpoints of the arc.
      if (xcv.test_flag(X_monotone_curve_2::IS_DIRECTED_RIGHT)) {
        // The arc is directed from left to right, so p becomes xcv1's target
        // and xcv2's source.
        xcv1.set_target(p);
        xcv2.set_source(p);

        if (! p.is_generating_conic(xcv.id())) {
          xcv1.target().set_generating_conic(xcv.id());
          xcv2.source().set_generating_conic(xcv.id());
        }
      }
      else {
        // The arc is directed from right to left, so p becomes xcv2's target
        // and xcv1's source.
        xcv1.set_source(p);
        xcv2.set_target(p);

        if (! p.is_generating_conic(xcv.id())) {
          xcv1.source().set_generating_conic(xcv.id());
          xcv2.target().set_generating_conic(xcv.id());
        }
      }
    }

  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object() const { return Split_2(*this); }

  class Intersect_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor.
     * \param traits The traits.
     */
    Intersect_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Find the intersections of the two given curves and insert them to the
     * given output iterator. As two segments may itersect only once, only a
     * single will be contained in the iterator.
     * \param cv1 The first arc.
     * \param cv2 The second arc.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv1,
                              const X_monotone_curve_2& xcv2,
                              OutputIterator oi) const
    { return intersect(xcv1, xcv2, m_traits.m_inter_map, oi); }

  private:
    /*! Compute the overlap with a given arc, which is supposed to have the same
     * supporting conic curve as this arc.
     * \param xcv1 The first arc.
     * \param xcv2 The second arc.
     * \param overlap Output: The overlapping arc (if any).
     * \return Whether we found an overlap.
     */
    bool compute_overlap(const X_monotone_curve_2& xcv1,
                         const X_monotone_curve_2& xcv2,
                         X_monotone_curve_2& overlap) const {
      // Check if the two arcs are identical.
      if (m_traits.equal_2_object()(xcv1,xcv2)) {
        overlap = xcv2;
        return true;
      }

      if (m_traits.is_strictly_between_endpoints(xcv1, xcv2.left())) {
        if (m_traits.is_strictly_between_endpoints(xcv1, xcv2.right())) {
          // Case 1 - *this:     +----------->
          //            arc:       +=====>
          overlap = xcv2;
          return true;
        }
        else {
          // Case 2 - *this:     +----------->
          //            arc:               +=====>
          overlap = xcv1;

          if (overlap.is_directed_right()) overlap.m_source = xcv2.left();
          else overlap.m_target = xcv2.left();

          return true;
        }
      }
      else if (m_traits.is_strictly_between_endpoints(xcv1, xcv2.right())) {
        // Case 3 - *this:     +----------->
        //            arc:   +=====>
        overlap = xcv1;
        if (overlap.is_directed_right()) overlap.m_target = xcv2.right();
        else overlap.m_source = xcv2.right();
        return true;
      }
      else if (m_traits.is_between_endpoints(xcv2, xcv1.source()) &&
               m_traits.is_between_endpoints(xcv2, xcv1.target()) &&
               (m_traits.is_strictly_between_endpoints(xcv2, xcv1.source()) ||
                m_traits.is_strictly_between_endpoints(xcv2, xcv1.target())))
      {
        // Case 4 - *this:     +----------->
        //            arc:   +================>
        overlap = xcv1;
        return true;
      }

      // If we reached here, there are no overlaps:
      return false;
    }

    /*! Intersect the supporing conic curves of this arc and the given arc.
     * \param arc The arc to intersect with.
     * \param inter_list The list of intersection points.
     */
    void intersect_supporting_conics(const X_monotone_curve_2& xcv1,
                                     const X_monotone_curve_2& xcv2,
                                     Intersection_list& inter_list) const {
      if (xcv1.is_special_segment() && ! xcv2.is_special_segment()) {
        // If one of the arcs is a special segment, make sure it is (arc).
        intersect_supporting_conics(xcv2, xcv1, inter_list);
        return;
      }

      const int deg1 =
        (xcv1.degree_mask() == X_monotone_curve_2::degree_1_mask()) ? 1 : 2;
      const int deg2 =
        (xcv2.degree_mask() == X_monotone_curve_2::degree_1_mask()) ? 1 : 2;
      const auto nt_traits = m_traits.m_nt_traits;
      Algebraic xs[4];
      int n_xs = 0;
      Algebraic ys[4];
      int n_ys = 0;

      if (xcv2.is_special_segment()) {
        // The second arc is a special segment (a*x + b*y + c = 0).
        if (xcv1.is_special_segment()) {
          // Both arc are sepcial segment, so they have at most one intersection
          // point.
          const auto* extra_data1 = xcv1.extra_data();
          const auto* extra_data2 = xcv2.extra_data();
          Algebraic denom =
            extra_data1->a * extra_data2->b - extra_data1->b * extra_data2->a;
          if (CGAL::sign (denom) != CGAL::ZERO) {
            xs[0] = (extra_data1->b * extra_data2->c -
                     extra_data1->c * extra_data2->b) / denom;
            n_xs = 1;

            ys[0] = (extra_data1->c * extra_data2->a -
                     extra_data1->a * extra_data2->c) / denom;
            n_ys = 1;
          }
        }
        else {
          const auto* extra_data2 = xcv2.extra_data();

          // Compute the x-coordinates of the intersection points.
          n_xs = compute_resultant_roots(*nt_traits,
                                         xcv1.alg_r(), xcv1.alg_s(),
                                         xcv1.alg_t(), xcv1.alg_u(),
                                         xcv1.alg_v(), xcv1.alg_w(),
                                         deg1,
                                         extra_data2->a,
                                         extra_data2->b,
                                         extra_data2->c,
                                         xs);
          CGAL_assertion(n_xs <= 2);

          // Compute the y-coordinates of the intersection points.
          n_ys = compute_resultant_roots(*nt_traits,
                                         xcv1.alg_s(), xcv1.alg_r(),
                                         xcv1.alg_t(), xcv1.alg_v(),
                                         xcv1.alg_u(), xcv1.alg_w(),
                                         deg1,
                                         extra_data2->b,
                                         extra_data2->a,
                                         extra_data2->c,
                                         ys);
          CGAL_assertion(n_ys <= 2);
        }
      }
      else {
        // Compute the x-coordinates of the intersection points.
        n_xs = compute_resultant_roots(*nt_traits,
                                       xcv1.r(), xcv1.s(), xcv1.t(),
                                       xcv1.u(), xcv1.v(), xcv1.w(),
                                       deg1,
                                       xcv2.r(), xcv2.s(), xcv2.t(),
                                       xcv2.u(), xcv2.v(), xcv2.w(),
                                       deg2,
                                       xs);
        CGAL_assertion(n_xs <= 4);

        // Compute the y-coordinates of the intersection points.
        n_ys = compute_resultant_roots(*nt_traits,
                                       xcv1.s(), xcv1.r(), xcv1.t(),
                                       xcv1.v(), xcv1.u(), xcv1.w(),
                                       deg1,
                                       xcv2.s(), xcv2.r(), xcv2.t(),
                                       xcv2.v(), xcv2.u(), xcv2.w(),
                                       deg2,
                                       ys);
        CGAL_assertion(n_ys <= 4);
      }

      // Pair the coordinates of the intersection points. As the vectors of
      // x and y-coordinates are sorted in ascending order, we output the
      // intersection points in lexicographically ascending order.
      Multiplicity mult;
      int i, j;

      if (xcv2.is_special_segment()) {
        if ((n_xs == 0) || (n_ys == 0)) return;

        if ((n_xs == 1) && (n_ys == 1)) {
          // Single intersection.
          Point_2 ip(xs[0], ys[0]);
          ip.set_generating_conic(xcv1.id());
          ip.set_generating_conic(xcv2.id());

          // In case the other curve is of degree 2, this is a tangency point.
          mult = ((deg1 == 1) || xcv1.is_special_segment()) ? 1 : 2;
          inter_list.push_back(Intersection_point(ip, mult));
        }
        else if ((n_xs == 1) && (n_ys == 2)) {
          Point_2 ip1(xs[0], ys[0]);
          ip1.set_generating_conic(xcv1.id());
          ip1.set_generating_conic(xcv2.id());

          inter_list.push_back(Intersection_point(ip1, 1));

          Point_2 ip2(xs[0], ys[1]);
          ip2.set_generating_conic(xcv1.id());
          ip2.set_generating_conic(xcv2.id());

          inter_list.push_back(Intersection_point(ip2, 1));
        }
        else if ((n_xs == 2) && (n_ys == 1)) {
          Point_2 ip1(xs[0], ys[0]);
          ip1.set_generating_conic(xcv1.id());
          ip1.set_generating_conic(xcv2.id());

          inter_list.push_back(Intersection_point(ip1, 1));

          Point_2 ip2(xs[1], ys[0]);
          ip2.set_generating_conic(xcv1.id());
          ip2.set_generating_conic(xcv2.id());

          inter_list.push_back(Intersection_point(ip2, 1));
        }
        else {
          CGAL_assertion((n_xs == 2) && (n_ys == 2));

          // The x-coordinates and the y-coordinates are given in ascending
          // order. If the slope of the segment is positive, we pair the
          // coordinates as is - otherwise, we swap the pairs.
          int ind_first_y(0), ind_second_y(1);

          const auto* extra_data2 = xcv2.extra_data();
          if (CGAL::sign(extra_data2->b) == CGAL::sign(extra_data2->a)) {
            ind_first_y = 1;
            ind_second_y = 0;
          }

          Point_2 ip1(xs[0], ys[ind_first_y]);
          ip1.set_generating_conic(xcv1.id());
          ip1.set_generating_conic(xcv2.id());

          inter_list.push_back(Intersection_point(ip1, 1));

          Point_2 ip2(xs[1], ys[ind_second_y]);
          ip2.set_generating_conic(xcv1.id());
          ip2.set_generating_conic(xcv2.id());

          inter_list.push_back(Intersection_point(ip2, 1));
        }

        return;
      }

      for (i = 0; i < n_xs; ++i) {
        for (j = 0; j < n_ys; ++j) {
          if (xcv1.is_on_supporting_conic(xs[i], ys[j]) &&
              xcv2.is_on_supporting_conic(xs[i], ys[j]))
          {
            // Create the intersection point and set its generating conics.
            Point_2 ip(xs[i], ys[j]);

            ip.set_generating_conic(xcv1.id());
            ip.set_generating_conic(xcv2.id());

            // Compute the multiplicity of the intersection point.
            if (deg1 == 1 && deg2 == 1) mult = 1;
            else mult = xcv1.multiplicity_of_intersection_point(xcv2, ip);

            // Insert the intersection point to the output list.
            inter_list.push_back(Intersection_point(ip, mult));
          }
        }
      }
    }

    /*! Compute the intersections with the given arc.
     * \param arc The given intersecting arc.
     * \param inter_map Maps conic pairs to lists of their intersection points.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator intersect(const X_monotone_curve_2& xcv1,
                             const X_monotone_curve_2& xcv2,
                             Intersection_map& inter_map,
                             OutputIterator oi) const
    {
      if (m_traits.has_same_supporting_conic(xcv1, xcv2)) {
        // Check for overlaps between the two arcs.
        X_monotone_curve_2 overlap;

        if (compute_overlap(xcv1, xcv2, overlap)) {
          // There can be just a single overlap between two x-monotone arcs:
          *oi++ = overlap;
          return oi;
        }

        // In case there is not overlap and the supporting conics are the same,
        // there cannot be any intersection points, unless the two arcs share
        // an end point.
        // Note that in this case we do not define the multiplicity of the
        // intersection points we report.
        auto alg_kernel = m_traits.m_alg_kernel;
        auto eq = alg_kernel->equal_2_object();
        if (eq(xcv1.left(), xcv2.left())) {
          Intersection_point ip(xcv1.left(), 0);
          *oi++ = ip;
        }

        if (eq(xcv1.right(), xcv2.right())) {
          Intersection_point ip(xcv1.right(), 0);
          *oi++ = ip;
        }

        if (eq(xcv1.left(), xcv2.right())) {
          Intersection_point ip(xcv1.left(), 0);
          *oi++ = ip;
        }

        if (eq(xcv1.right(), xcv2.left())) {
          Intersection_point ip(xcv1.right(), 0);
          *oi++ = ip;
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

      if (xcv1.id().is_valid() && xcv2.id().is_valid()) {
        if (xcv1.id() < xcv2.id()) conic_pair = Conic_pair(xcv1.id(), xcv2.id());
        else conic_pair = Conic_pair(xcv2.id(), xcv1.id());
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
        intersect_supporting_conics(xcv1, xcv2, inter_list);
        if (! invalid_ids) inter_map[conic_pair] = inter_list;
      }
      else {
        // Obtain the precomputed intersection points from the map.
        inter_list = (*map_iter).second;
      }

      // Go over the list of intersection points and report those that lie on
      // both x-monotone arcs.
      for (auto iter = inter_list.begin(); iter != inter_list.end(); ++iter) {
        if (m_traits.is_between_endpoints(xcv1, (*iter).first) &&
            m_traits.is_between_endpoints(xcv2, (*iter).first))
        {
          *oi++ = *iter;
        }
      }

      return oi;
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const { return Intersect_2(*this); }

  class Are_mergeable_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits (in case it has state)
     */
    Are_mergeable_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param xcv1 The first arc.
     * \param xcv2 The second arc.
     * \return `true` if the two curves are mergeable; that is, they are
     *         supported by the same curve and share a common endpoint);
     *         `false` otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    { return can_merge_with(xcv1, xcv2); }

  private:
    /*! Check whether it is possible to merge the arc with the given arc.
     * \param xcv1 The first arc.
     * \param xcv2 The second arc.
     * \return `true` if it is possible to merge the two arcs;
     *         `false` otherwise.
     */
    bool can_merge_with(const X_monotone_curve_2& xcv1,
                        const X_monotone_curve_2& xcv2) const {
      // In order to merge the two arcs, they should have the same supporting
      // conic.
      if (! m_traits.has_same_supporting_conic(xcv1, xcv2)) return false;

      // Check if the left endpoint of one curve is the right endpoint of the
      // other.
      auto eq = m_traits.m_alg_kernel->equal_2_object();
      return (eq(xcv1.right(), xcv2.left()) || eq(xcv1.left(), xcv2.right()));
    }
  };

  /*! Obtain an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(*this); }

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits (in case it has state)
     */
    Merge_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Merge two given x-monotone curves into a single curve (segment).
     * \param xcv1 The first arc.
     * \param xcv2 The second arc.
     * \param xcv The merged arc.
     * \pre The two arcs are mergeable.
     */
    void operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2,
                    X_monotone_curve_2& xcv) const
    {
      CGAL_precondition(m_traits.are_mergeable_2_object()(xcv2, xcv1));
      xcv = xcv1;
      merge(xcv, xcv2);
    }

  private:
    /*! Merge the current arc with the given arc.
     * \param xcv1 The first arc to merge with.
     * \param xcv2 The second arc to merge with.
     * \pre The two arcs are mergeable.
     */
    void merge(X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2) const {
      // Check whether we should extend the arc to the left or to the right.
      auto eq = m_traits.m_alg_kernel->equal_2_object();
      if (eq(xcv1.right(), xcv2.left())) {
        // Extend the arc to the right.
        if (xcv1.test_flag(X_monotone_curve_2::IS_DIRECTED_RIGHT))
          xcv1.set_target(xcv2.right());
        else xcv1.set_source(xcv2.right());
      }
      else {
        CGAL_precondition(eq(xcv1.left(), xcv2.right()));

        // Extend the arc to the left.
        if (xcv1.test_flag(X_monotone_curve_2::IS_DIRECTED_RIGHT))
          xcv1.set_source(xcv2.left());
        else xcv1.set_target(xcv2.left());
      }
    }
  };

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object() const { return Merge_2(*this); }

  //@}

  /*! \name Auxiliary Functor definitions, used gor, e.g., the landmarks
   * point-location strategy and the drawing function.
   */
  //@{
  typedef double                                        Approximate_number_type;
  typedef CGAL::Cartesian<Approximate_number_type>      Approximate_kernel;
  typedef Approximate_kernel::Point_2                   Approximate_point_2;

  class Approximate_curve_length_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits.
     */
    Approximate_curve_length_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Obtain an approximation of the length of a curve.
     * \param xcv The curve.
     * \return An approximation of the length of `xcv`.
     */
    Approximate_number_type operator()(const X_monotone_curve_2& xcv) const {
      if (xcv.orientation() == COLLINEAR) return segment_length(xcv);
      CGAL::Sign sign_conic = CGAL::sign(4*xcv.r()*xcv.s() - xcv.t()*xcv.t());
      if (sign_conic == POSITIVE) return ellipse_length(xcv);
      if (sign_conic == NEGATIVE) return hyperbola_length(xcv);
      return parabola_length(xcv);
    }

  private:
    /*! Obtain the segment length.
     */
    double segment_length(const X_monotone_curve_2& xcv) {
      auto min_vertex = m_traits.construct_min_vertex_2_object();
      auto max_vertex = m_traits.construct_max_vertex_2_object();
      const auto& minv = min_vertex(xcv);
      const auto& maxv = max_vertex(xcv);
      auto x1 = CGAL::to_double(minv.x());
      auto y1 = CGAL::to_double(minv.y());
      auto x2 = CGAL::to_double(maxv.x());
      auto y2 = CGAL::to_double(maxv.y());
      auto dx = x2 - x1;
      auto dy = y2 - y1;
      double l = std::sqrt(x2*x2 + y2*y2);
      return l;
    }

    /*! The formula for the arc length of a parabola is:
     * L = 1/2 * √(b^2+16⋅a^2) + b^2/(8*a) * ln((4*a+√(b^2+16⋅a^2))/b)
     * where:
     * L is the length of the parabola arc.
     * a is the length along the parabola axis.
     * b is the length of the chord perpendicular to the axis.
     *
     *      ---
     *     / | \
     *    /  |a \
     *   /   |   \
     *  /---------\
     * /    b      \
     *
     */
    double parabolic_arc_length(double a, double b) const {
      if (a == 0) return b;
      if (b == 0) return a;
      auto b_sqr = b*b;
      auto tmp = std::sqrt(b_sqr+16.0*a*a);
      return tmp/2.0 + b_sqr*std::log((4.0*a + tmp)/b)/(8.0*a);
    }

    /*! Obtain the parabolic arc length.
     */
    double parabola_length(const X_monotone_curve_2& xcv) {
      double r_m, t_m, s_m, u_m, v_m, w_m;
      double cost, sint;
      double xs_t, ys_t, xt_t, yt_t;
      double a;
      double ts, tt;
      double cx, cy;
      m_traits.approximate_parabola(xcv,
                                    r_m, t_m, s_m, u_m, v_m, w_m, cost, sint,
                                    xs_t, ys_t, xt_t, yt_t, a, ts, tt, cx, cy);

      auto ds = parabolic_arc_length(xs_t, 2.0*std::abs(ys_t));
      auto dt = parabolic_arc_length(xt_t, 2.0*std::abs(yt_t));
      auto d = (CGAL::sign(ys_t) == CGAL::sign(yt_t)) ?
        std::abs(ds - dt)/2.0 : (ds + dt)/2.0;
      // std::cout << "d, ds, dt = " << d << ", " << ds << "," << dt
      //           << std::endl;
      return d;
    }

    double ellipse_length(const X_monotone_curve_2& xcv) {
      double r_m, t_m, s_m, u_m, v_m, w_m;
      double cost, sint;
      double xs_t, ys_t, xt_t, yt_t;
      double a, b;
      double cx, cy;
      double ts, tt;
      m_traits.approximate_ellipse(xcv, r_m, t_m, s_m, u_m, v_m, w_m, cost, sint,
                                   xs_t, ys_t, ts, xt_t, yt_t, tt,
                                   a, b, cx, cy);

      namespace bm = boost::math;
      auto ratio = b/a;
      auto k = std::sqrt(1 - (ratio*ratio));
      auto ds = a*bm::ellint_2(k, ts);
      auto dt = a*bm::ellint_2(k, tt);
      auto d = std::abs(dt - ds);
      // std::cout << "d,ds,dt: " << d << "," << ds << ", " << dt << std::endl;
      return d;
    }

    double hyperbola_length(const X_monotone_curve_2& /* xcv */) {
      CGAL_error_msg("Not implemented yet!");
      double l(0.0);
      return l;
    }
  };

  class Approximate_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits.
     */
    Approximate_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Obtain an approximation of a point coordinate.
     * \param p The exact point.
     * \param i The coordinate index (either 0 or 1).
     * \pre i is either 0 or 1.
     * \return An approximation of p's x-coordinate (if i == 0), or an
     *         approximation of p's y-coordinate (if i == 1).
     */
    Approximate_number_type operator()(const Point_2& p, int i) const {
      CGAL_precondition((i == 0) || (i == 1));

      if (i == 0) return CGAL::to_double(p.x());
      else return CGAL::to_double(p.y());
    }

    /*! Obtain an approximation of a point.
     */
    Approximate_point_2 operator()(const Point_2& p) const
    { return Approximate_point_2(operator()(p, 0), operator()(p, 1)); }

    /*! Obtain an approximation of an \f$x\f$-monotone curve.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv, double error,
                              OutputIterator oi, bool l2r = true) const {
      if (xcv.orientation() == COLLINEAR)
        return approximate_segment(xcv, oi, l2r);
      CGAL::Sign sign_conic = CGAL::sign(4*xcv.r()*xcv.s() - xcv.t()*xcv.t());
      if (sign_conic == POSITIVE)
        return approximate_ellipse(xcv, error, oi, l2r);
      if (sign_conic == NEGATIVE)
        return approximate_hyperbola(xcv, error, oi, l2r);
      return approximate_parabola(xcv, error, oi, l2r);
    }

  private:
    /*! Handle segments.
     */
    template <typename OutputIterator>
    OutputIterator approximate_segment(const X_monotone_curve_2& xcv,
                                       OutputIterator oi, bool l2r) const {
      // std::cout << "SEGMENT\n";
      auto min_vertex = m_traits.construct_min_vertex_2_object();
      auto max_vertex = m_traits.construct_max_vertex_2_object();
      const auto& src = (l2r) ? min_vertex(xcv) : max_vertex(xcv);
      const auto& trg = (l2r) ? max_vertex(xcv) : min_vertex(xcv);
      auto xs = CGAL::to_double(src.x());
      auto ys = CGAL::to_double(src.y());
      auto xt = CGAL::to_double(trg.x());
      auto yt = CGAL::to_double(trg.y());
      *oi++ = Approximate_point_2(xs, ys);
      *oi++ = Approximate_point_2(xt, yt);
      return oi;
    }

    /*! Transform a point. In particular, rotate the canonical point
     * (`xc`,`yc`) by an angle, the sine and cosine of which are `sint` and
     * `cost`, respectively, and translate by (`cx`,`cy`).
     */
    void transform_point(double xc, double yc,
                         double cost, double sint,
                         double cx, double cy,
                         double& x, double& y) const {
      x = xc*cost - yc*sint + cx;
      y = xc*sint + yc*cost + cy;
    }

    /*! Handle ellipses.
     * The general equation of an ellipse is:
     *   r·𝑥^2 + s·𝑦^2 + t·𝑥·𝑦 + u·𝑥 + v·𝑦 + w = 0
     * where 4·r·s−t^2 > 0
     * We eliminate t so that the x·y term vanishes, applying an inverse
     * rotation. Then, we compute the radi and the center. Finaly, we rotate
     * back. The angle of rotation is given by:
     *   𝑡𝑎𝑛(2𝜃) = 𝐵 / (𝐴−𝐶)
     * Then
     *   𝑐𝑜𝑠(2𝜃) = sqrt((1 + 𝑐𝑜𝑠(2𝜃))/2)
     *   𝑠𝑖𝑛(2𝜃) = sqrt((1 - 𝑐𝑜𝑠(2𝜃))/2)
     * The coefficients of the new ellipse are given by:
     *   r′ = r·𝑐𝑜𝑠(𝜃)^2 + s·𝑠𝑖𝑛(𝜃)^2 + t·𝑐𝑜𝑠(𝜃)·𝑠𝑖𝑛(𝜃)
     *   s′ = r·𝑠𝑖𝑛(𝜃)^2 + s·𝑐𝑜𝑠(𝜃)^2 − t·𝑐𝑜𝑠(𝜃)·𝑠𝑖𝑛(𝜃)
     *   t′ = 0
     *   u′ = u·𝑐𝑜𝑠(𝜃) + v·𝑠𝑖𝑛(𝜃)
     *   v′ = −u·𝑠𝑖𝑛(𝜃) + v·𝑐𝑜𝑠(𝜃)
     *   w′ = w
     *  After writing this equation in the form:
     *   (𝑥′ − 𝐶𝑥′)^2    (𝑦′ − 𝐶𝑦′)^2
     *   ----------- + ------------ = 1
     *       𝑎^2            𝑏^2
     * We get:
     *          -u′
     *   𝐶𝑥′ = ----
     *         2·r′
     *
     *          -v′
     *   𝐶𝑦′ = ----
     *         2·s′
     *
     *         −4·w′·r′·s′ + s′·u′2 + r′·v′2
     *   a^2 = ----------------------------
     *          4·r′^2·s′
     *
     *         −4·w′·r′·s′ + s′·u′2 + r′·v′2
     *   b^2 = ----------------------------
     *         4·r′·s′^2
     *
     * Rotate back about angle 𝜃 to find the coordinates of the center:
     *   𝐶𝑥 = 𝐶𝑥′·𝑐𝑜𝑠𝜃 − 𝐶𝑦′·𝑠𝑖𝑛𝜃
     *   𝐶𝑦 = 𝐶𝑥′·𝑠𝑖𝑛𝜃 + 𝐶𝑦′·𝑐𝑜𝑠𝜃
     *
     * The parametric formula of an ellipse centered at the origin with major
     * axis parallel to the x-axis and minor axis parallel to the y-axis is:
     *  𝑥(𝛼) = a·𝑐𝑜𝑠(𝛼)
     *  𝑦(𝛼) = b·𝑠𝑖𝑛(𝛼)
     * where a is the major radius and b is the minor radius.
     *
     * The rotation transformation is fiven by
     *  𝑥 = 𝑥(𝑡)·𝑐𝑜𝑠(𝜃) − 𝑦(𝑡)·𝑠𝑖𝑛(𝜃)
     *  𝑦 = 𝑥(𝑡)·𝑠𝑖𝑛(𝜃) + 𝑦(𝑡)·𝑐𝑜𝑠(𝜃)
     * Where 𝜃 is the rotation angle
     *
     * Combining the above we get
     *  𝑥(𝛼) = a·𝑐𝑜𝑠(𝛼)·𝑐𝑜𝑠(𝜃) − b·𝑠𝑖𝑛(𝛼)·𝑠𝑖𝑛(𝜃)
     *  𝑦(𝛼) = a·𝑐𝑜𝑠(𝛼)·𝑠𝑖𝑛(𝜃) + b·𝑠𝑖𝑛(𝛼)·𝑐𝑜𝑠(𝜃)
     *
     * To shift from the center we add 𝐶𝑥 to the 𝑥(𝛼) and 𝐶𝑦 to 𝑦(𝛼).
     * Therefore, the equations of a Rotated Ellipse are:
     *   𝑥(𝛼) = a·𝑐𝑜𝑠(𝛼)·𝑐𝑜𝑠(𝜃) − b·𝑠𝑖𝑛(𝛼)·𝑠𝑖𝑛(𝜃) + 𝐶𝑥
     *   𝑦(𝛼) = a·𝑐𝑜𝑠(𝛼)·𝑠𝑖𝑛(𝜃) + b·𝑠𝑖𝑛(𝛼)·𝑐𝑜𝑠(𝜃) + 𝐶𝑦
     *
     * @param error the error bound of the generated approximation. This is
     *              the Hausdorff distance between the arc and the polyline,
     *              which approximates the arc.
     */
    template <typename OutputIterator>
    OutputIterator approximate_ellipse(const X_monotone_curve_2& xcv,
                                       double error, OutputIterator oi,
                                       bool l2r = true) const {
      // std::cout << "ELLIPSE\n";
      auto min_vertex = m_traits.construct_min_vertex_2_object();
      auto max_vertex = m_traits.construct_max_vertex_2_object();
      const auto& src = (l2r) ? min_vertex(xcv) : max_vertex(xcv);
      const auto& trg = (l2r) ? max_vertex(xcv) : min_vertex(xcv);
      auto xs = CGAL::to_double(src.x());
      auto ys = CGAL::to_double(src.y());
      auto xt = CGAL::to_double(trg.x());
      auto yt = CGAL::to_double(trg.y());
      // std::cout << "curve: (" << xs << "," << ys
      //           << ") => (" << xt << "," << yt << ")"
      //           << std::endl;

      double r_m, t_m, s_m, u_m, v_m, w_m;
      double cost, sint;
      double xs_t, ys_t, xt_t, yt_t;
      double a, b;
      double cx, cy;
      double ts, tt;
      m_traits.approximate_ellipse(xcv, r_m, t_m, s_m, u_m, v_m, w_m, cost, sint,
                                   xs_t, ys_t, ts, xt_t, yt_t, tt,
                                   a, b, cx, cy, l2r);
      // std::cout << "a, b: " << a << "," << b << std::endl;

      *oi++ = Approximate_point_2(xs, ys);
      add_points(xs_t, ys_t, ts, xt_t, yt_t, tt, error, oi,
                 [&](double tm, double& xm, double& ym) {
                   elliptic_point(a, b, tm, xm, ym);
                 },
                 [&](double xc, double& yc, double& x, double& y) {
                   transform_point(xc, yc, cost, sint, cx, cy, x, y);
                 });
      *oi++ = Approximate_point_2(xt, yt);
      return oi;
    }

    /*! Add either an elliptic or a hyperbilc point.
     * The arc endpoints are (`x1`, `y1`) and (`x2`, `y2`).
     * In our parametric representations for ellipses and hyperbolas the
     * following holds:
     *   p1 = (x1, y1); x1 = x(t1), y1 = y(t1), and
     *   p2 = (x2, y2); x2 = x(t2), y2 = y(t2)
     * The Hausdorff distance between the arc and the segment (p1,p2) is
     * at (xm,ym), where xm = x(tm), ym = y(tmp), and tm = (t1 + t2) / 2.
     * \param[in] x1 the canonical-arc source-point \f$x\f$-coordinate.
     * \param[in] y1 the canonical-arc source-point \f$y\f$-coordinate.
     * \param[in] t1 the source-point parameter; \f$yx1 = x(t1) and y1 = y(t1)\f$y.
     * \param[in] x2 the canonical-arc target-point \f$x\f$-coordinate.
     * \param[in] y2 the canonical arc target-point \f$y\f$-coordinate.
     * \param[in] t2 the target-point parameter; \f$yx2 = x(t2) and y2 = y(t2)\f$y.
     * \param[in] error
     * \param[out] oi
     * \param[in] op a function that computes a point \f$(x(t),y(t))\f$ given \f$t\f$.
     * \param[in] transform a function that transforms a canonical point to an
     *                      actual point
     *
     * Observe that in our parametric representation for parabolas, the
     * expression for tm is different.
     */
    template <typename OutputIterator, typename Op, typename Transform>
    OutputIterator add_points(double x1, double y1, double t1,
                              double x2, double y2, double t2,
                              double error, OutputIterator oi,
                              Op op, Transform transform) const {
      auto tm = (t1 + t2)*0.5;

      // Compute the canocal point where the error is maximal.
      double xm, ym;
      op(tm, xm, ym);

      auto dx = x2 - x1;
      auto dy = y2 - y1;

      // Compute the error; abort if it is below the threshold
      auto l = std::sqrt(dx*dx + dy*dy);
      auto e = std::abs((xm*dy - ym*dx + x2*y1 - x1*y2) / l);
      if (e < error) return oi;

      double x, y;
      transform(xm, ym, x, y);
      add_points(x1, y1, t1, xm, ym, tm, error, oi, op, transform);
      *oi++ = Approximate_point_2(x, y);
      add_points(xm, ym, tm, x2, y2, t2, error, oi, op, transform);
      return oi;
    }

    /*! Compute the elliptic point given the parameter t and the transform
     * data, that is, the center (translation) and the sin and cos of the
     * rotation angle.
     */
    void elliptic_point(double a, double b, double t,
                        double& x, double& y) const {
      x = a * std::cos(t);
      y = b * std::sin(t);
    }

    /*! Handle parabolas.
     * The arc-length closed form can be found here:
     * https://www.vcalc.com/wiki/vCalc/Parabola+-+arc+length
     */
    template <typename OutputIterator>
    OutputIterator approximate_parabola(const X_monotone_curve_2& xcv,
                                        double error, OutputIterator oi,
                                        bool l2r = true)
      const {
      // std::cout << "PARABOLA\n";
      auto min_vertex = m_traits.construct_min_vertex_2_object();
      auto max_vertex = m_traits.construct_max_vertex_2_object();
      const auto& src = (l2r) ? min_vertex(xcv) : max_vertex(xcv);
      const auto& trg = (l2r) ? max_vertex(xcv) : min_vertex(xcv);
      auto xs = CGAL::to_double(src.x());
      auto ys = CGAL::to_double(src.y());
      auto xt = CGAL::to_double(trg.x());
      auto yt = CGAL::to_double(trg.y());
      // std::cout << "curve: (" << xs << "," << ys
      //           << ") => (" << xt << "," << yt << ")"
      //           << std::endl;

      double r_m, t_m, s_m, u_m, v_m, w_m;
      double cost, sint;
      double xs_t, ys_t, xt_t, yt_t;
      double a;
      double ts, tt;
      double cx, cy;
      m_traits.approximate_parabola(xcv,
                                    r_m, t_m, s_m, u_m, v_m, w_m, cost, sint,
                                    xs_t, ys_t, xt_t, yt_t, a, ts, tt, cx, cy,
                                    l2r);
      // std::cout << "sint, cost: " << sint << "," << cost << std::endl;
      // std::cout << "a: " << a << std::endl;
      // std::cout << "xs' = " << xs_t << "," << ys_t << std::endl;
      // std::cout << "xt' = " << xt_t << "," << yt_t << std::endl;
      // std::cout << "ts,tt = " << ts << "," << tt << std::endl;

      *oi++ = Approximate_point_2(xs, ys);
      add_parabolic_points(xs_t, ys_t, ts, xt_t, yt_t, tt, error, oi,
                           [&](double tm, double& xm, double& ym) {
                             parabolic_point(a, tm, xm, ym);
                           },
                           [&](double xc, double& yc, double& x, double& y) {
                             transform_point(xc, yc, cost, sint, cx, cy, x, y);
                           });
      *oi++ = Approximate_point_2(xt, yt);
      return oi;
    }

    /*! Add either an elliptic or a hyperbilc point.
     * The arc endpoints are (`x1`, `y1`) and (`x2`, `y2`).
     * In our parametric representations for ellipses and hyperbolas the
     * following holds:
     *   \f$p1 = (x1, y1); x1 = x(t1), y1 = y(t1)\f$, and
     *   \f$p2 = (x2, y2); x2 = x(t2), y2 = y(t2)\f$
     * The Hausdorff distance between the arc and the segment (p1,p2) is
     * at (xm,ym), where xm = x(tm), ym = y(tmp), and tm = (x2-x1) / (y2-y1).
     * \param[in] x1 the canonical-arc source-point \f$x\f$-coordinate.
     * \param[in] y1 the canonical-arc source-point \f$y\f$-coordinate.
     * \param[in] t1 the source-point parameter; \f$yx1 = x(t1) and y1 = y(t1)\f$y.
     * \param[in] x2 the canonical-arc target-point \f$x\f$-coordinate.
     * \param[in] y2 the canonical arc target-point \f$y\f$-coordinate.
     * \param[in] t2 the target-point parameter; \f$yx2 = x(t2) and y2 = y(t2)\f$y.
     * \param[in] error
     * \param[out] oi
     * \param[in] op a function that computes a point \f$(x(t),y(t))\f$ given \f$t\f$.
     * \param[in] transform a function that transforms a canonical point to an
     *                      actual point
     *
     * Observe that in our parametric representation for ellipses and
     * hyperbolas, the expression for tm is different.
     */
    template <typename OutputIterator, typename Op, typename Transform>
    OutputIterator add_parabolic_points(double x1, double y1, double t1,
                                        double x2, double y2, double t2,
                                        double error, OutputIterator oi,
                                        Op op, Transform transform) const {
      auto dx = x2 - x1;
      auto dy = y2 - y1;
      auto tm = (dy == 0) ? 0 : dx / dy;

      // Compute the canocal point where the error is maximal.
      double xm, ym;
      op(tm, xm, ym);

      // Compute the error and abort if it is below the threshold
      auto l = std::sqrt(dx*dx + dy*dy);
      auto e = std::abs((xm*dy - ym*dx + x2*y1 - x1*y2) / l);
      if (e < error) return oi;

      // Compute the actual (transformed) point
      double x, y;
      transform(xm, ym, x, y);
      add_parabolic_points(x1, y1, t1, xm, ym, tm, error, oi, op, transform);
      *oi++ = Approximate_point_2(x, y);
      add_parabolic_points(xm, ym, tm, x2, y2, t2, error, oi, op, transform);
      return oi;
    }

    /*! Compute the parabolic point given the parameter t and the transform
     * data, that is, the center (translation) and the sin and cos of the
     * rotation angle.
     */
    void parabolic_point(double a, double t, double& x, double& y) const {
      x = a*t*t;
      y = 2.0*a*t;
    }

    /*! Handle hyperbolas.
     */
    template <typename OutputIterator>
    OutputIterator approximate_hyperbola(const X_monotone_curve_2& xcv,
                                         double error, OutputIterator oi,
                                         bool l2r = true) const {
      // std::cout << "HYPERBOLA\n";
      auto min_vertex = m_traits.construct_min_vertex_2_object();
      auto max_vertex = m_traits.construct_max_vertex_2_object();
      const auto& src = (l2r) ? min_vertex(xcv) : max_vertex(xcv);
      const auto& trg = (l2r) ? max_vertex(xcv) : min_vertex(xcv);
      auto xs = CGAL::to_double(src.x());
      auto ys = CGAL::to_double(src.y());
      auto xt = CGAL::to_double(trg.x());
      auto yt = CGAL::to_double(trg.y());

      double r_m, t_m, s_m, u_m, v_m, w_m;
      double cost, sint;
      double xs_t, ys_t, xt_t, yt_t;
      double a, b;
      double cx, cy;
      double ts, tt;
      m_traits.approximate_hyperbola(xcv, r_m, t_m, s_m, u_m, v_m, w_m,
                                     cost, sint,
                                     xs_t, ys_t, ts, xt_t, yt_t, tt,
                                     a, b, cx, cy, l2r);
      // std::cout << "a, b: " << a << "," << b << std::endl;
      // std::cout << "ts, tt: " << ts << "," << tt << std::endl;

      // std::cout << "a, b: " << a << "," << b << std::endl;
      *oi++ = Approximate_point_2(xs, ys);
      add_points(xs_t, ys_t, ts, xt_t, yt_t, tt, error, oi,
                 [&](double tm, double& xm, double& ym) {
                   hyperbolic_point(a, b, tm, xm, ym);
                 },
                 [&](double xc, double& yc, double& x, double& y) {
                   transform_point(xc, yc, cost, sint, cx, cy, x, y);
                 });
      *oi++ = Approximate_point_2(xt, yt);
      return oi;
    }

    /*! Compute the hyperbolic point given the parameter t and the transform
     * data, that is, the center (translation) and the sin and cos of the
     * rotation angle.
     */
    void hyperbolic_point(double a, double b, double t,
                          double& x, double& y) const {
      x = a * std::cosh(t);
      y = b * std::sinh(t);
    }
  };

  /*! Obtain an Approximate_2 functor object. */
  Approximate_2 approximate_2_object() const { return Approximate_2(*this); }

  //! Functor
  class Construct_x_monotone_curve_2 {
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits.
     */
    Construct_x_monotone_curve_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Construct an \f$x\f$-monotone arc from a conic arc.
     * \param cv The given curve.
     * \pre cv is \f$x\f$-monotone.
     */
    X_monotone_curve_2 operator()(const Curve_2& cv) const {
      CGAL_precondition(cv.is_valid() && is_x_monotone(cv));
      X_monotone_curve_2 xcv(cv);
      m_traits.set_x_monotone(xcv);
      return xcv;
    }

    /*! Construct an \f$x\f$-monotone arc from a conic arc.
     * \param xcv The given curve.
     * \param id The ID of the base curve.
     */
    X_monotone_curve_2 operator()(const Curve_2& cv, const Conic_id& id) const {
      X_monotone_curve_2 xcv(cv, id);
      CGAL_precondition(xcv.is_valid() && id.is_valid());
      m_traits.set_x_monotone(xcv);
      return xcv;
    }

    /*! Construct an \f$x\f$-monotone sub-arc from a conic arc.
     * \param cv The given (base) arc.
     * \param source The source point.
     * \param target The target point.
     * \param id The id of the base arc.
     */
    X_monotone_curve_2 operator()(const Curve_2& cv,
                                  const Point_2& source, const Point_2& target,
                                  const Conic_id& id) const
    {
      // Set the two endpoints.
      X_monotone_curve_2 xcv(cv, id);
      xcv.set_source(source);
      xcv.set_target(target);
      CGAL_precondition(xcv.is_valid() && id.is_valid());
      m_traits.set_x_monotone(xcv);
      return xcv;
    }

    /*! Return an \f$x\f$-monotone curve connecting the two given endpoints.
     * \param source The first point.
     * \param target The second point.
     * \pre `source` and `target` must not be the same.
     * \return A segment connecting `source` and `target`.
     */
    X_monotone_curve_2 operator()(const Point_2& source, const Point_2& target)
      const
    {
      X_monotone_curve_2 xcv;

      // Set the basic properties.
      xcv.set_endpoints(source, target);
      xcv.set_orientation(COLLINEAR);
      xcv.set_flag(Curve_2::IS_VALID);

      // Set the other properties.
      xcv.set_flag(X_monotone_curve_2::DEGREE_1);
      xcv.set_flag(X_monotone_curve_2::IS_SPECIAL_SEGMENT);
      xcv.update_extra_data();

      auto cmp_xy = m_traits.m_alg_kernel->compare_xy_2_object();
      Comparison_result res = cmp_xy(source, target);
      CGAL_precondition(res != EQUAL);
      if (res == SMALLER) xcv.set_flag(X_monotone_curve_2::IS_DIRECTED_RIGHT);

      // Check whether the segment is vertical.
      if (CGAL::sign(xcv.extra_data()->b) == ZERO)
        xcv.set_flag(X_monotone_curve_2::IS_VERTICAL_SEGMENT);

      return xcv;
    }

    /*! Construct a special segment of a given line connecting to given
     * endpoints.
     * \param a, b, c The coefficients of the supporting line (`ax + by + c = 0`).
     * \param source The source point.
     * \param target The target point.
     */
    X_monotone_curve_2 operator()(const Algebraic& a, const Algebraic& b,
                                  const Algebraic& c,
                                  const Point_2& source, const Point_2& target)
      const
    {
      auto cmp_xy = m_traits.m_alg_kernel->compare_xy_2_object();
      Comparison_result res = cmp_xy(source, target);
      CGAL_precondition(res != EQUAL);

      X_monotone_curve_2 xcv;
      // Make sure the two endpoints lie on the supporting line.
      CGAL_precondition(CGAL::sign(a*source.x()+b*source.y()+c) == CGAL::ZERO);
      CGAL_precondition(CGAL::sign(a*target.x()+b*target.y()+c) == CGAL::ZERO);

      // Set the basic properties.
      xcv.set_endpoints(source, target);
      xcv.set_orientation(COLLINEAR);
      xcv.set_flag(Curve_2::IS_VALID);

      // Set the other properties.
      if (res == SMALLER) xcv.set_flag(X_monotone_curve_2::IS_DIRECTED_RIGHT);

      xcv.set_flag(X_monotone_curve_2::DEGREE_1);
      xcv.set_flag(X_monotone_curve_2::IS_SPECIAL_SEGMENT);
      xcv.set_extra_data(a, b, c, ZERO);

      // Check whether the segment is vertical.
      if (CGAL::sign(xcv.extra_data()->b) == ZERO)
        xcv.set_flag(X_monotone_curve_2::IS_VERTICAL_SEGMENT);

      return xcv;
    }

  private:
    /*! Determine whether the arc is \f$x\f$-monotone.
     */
    bool is_x_monotone(const Curve_2& cv) const {
      // Collect vertical tangency points.
      Alg_point_2 vtan_ps[2];
      auto res = m_traits.vertical_tangency_points(cv, vtan_ps);
      return (res == 0);
    }

    /*! Determine whether the arc is \f$y\f$-monotone.
     */
    bool is_y_monotone(const Curve_2& cv) const {
      // Collect horizontal tangency points.
      Alg_point_2 htan_ps[2];
      auto res = m_traits.horizontal_tangency_points(cv, htan_ps);
      return (res == 0);
    }

  };

  /*! Obtain a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object () const
  { return Construct_x_monotone_curve_2(*this); }

  //! Constructor of conic arcs
  class Construct_curve_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits.
     */
    Construct_curve_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Construct an empty curve.
     */
    Curve_2 operator()() const { return Curve_2(); }

    /*! Construct a conic arc which is the full conic:
     *   `C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0`
     * \pre The conic C must be an ellipse (so `4rs - t^2 > 0`).
     */
    Curve_2 operator()(const Rational& r, const Rational& s, const Rational& t,
                       const Rational& u, const Rational& v, const Rational& w)
      const
    {
      // Ensure that the given curve is an ellipse (4rs - t^2 is positive).
      CGAL_precondition(CGAL::sign(4*r*s - t*t) == POSITIVE);

      // Set the arc to be the full conic (and compute the orientation).
      Rational rat_coeffs[6];
      rat_coeffs[0] = r;
      rat_coeffs[1] = s;
      rat_coeffs[2] = t;
      rat_coeffs[3] = u;
      rat_coeffs[4] = v;
      rat_coeffs[5] = w;
      Curve_2 arc;
      m_traits.set_full(arc, rat_coeffs, true);
      return arc;
    }

    /*! Construct a conic arc that lies on the conic:
     *   `C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0`
     * \param orient The orientation of the arc (clockwise or counterclockwise).
     * \param source The source point.
     * \param target The target point.
     * \pre The source and the target must be on the conic boundary and must
     * not be the same.
     */
    Curve_2 operator()(const Rational& r, const Rational& s, const Rational& t,
                       const Rational& u, const Rational& v, const Rational& w,
                       Orientation orient,
                       const Point_2& source, const Point_2& target) const {
      // Make sure that the source and the taget are not the same.
      const auto alg_kernel = m_traits.m_alg_kernel;
      CGAL_precondition_code(auto eq = alg_kernel->compare_xy_2_object());
      CGAL_precondition(eq(source, target) != EQUAL);
      // Set the arc properties (no need to compute the orientation).
      Rational rat_coeffs[6] = {r, s, t, u, v, w};
      Curve_2 arc;
      arc.set_orientation(orient);
      arc.set_endpoints(source, target);
      m_traits.set(arc, rat_coeffs);
      return arc;
    }

    /*! Construct a conic arc that is a circular arc from given three points.
     * \param p1 The arc source.
     * \param p2 A point in the interior of the arc.
     * \param p3 The arc target.
     * \pre The three points must not be collinear.
     */
    Curve_2 operator()(const Rat_point_2& p1, const Rat_point_2& p2,
                       const Rat_point_2& p3) const {
      Curve_2 arc;

      // Set the source and target.
      const Rational& x1 = p1.x();
      const Rational& y1 = p1.y();
      const Rational& x2 = p2.x();
      const Rational& y2 = p2.y();
      const Rational& x3 = p3.x();
      const Rational& y3 = p3.y();

      const auto nt_traits = m_traits.m_nt_traits;
      const auto alg_kernel = m_traits.m_alg_kernel;
      Point_2 source(nt_traits->convert(x1), nt_traits->convert(y1));
      Point_2 target(nt_traits->convert(x3), nt_traits->convert(y3));
      arc.set_endpoints(source, target);

      // Make sure that the source and the taget are not the same.
      CGAL_precondition_code(auto cmp_xy = alg_kernel->compare_xy_2_object());
      CGAL_precondition(cmp_xy(source, target) != EQUAL);

      // Compute the lines: A1*x + B1*y + C1 = 0,
      //               and: A2*x + B2*y + C2 = 0,
      // where:
      const Rational two(2);

      const Rational A1 = two*(x1 - x2);
      const Rational B1 = two*(y1 - y2);
      const Rational C1 = y2*y2 - y1*y1 + x2*x2 - x1*x1;

      const Rational A2 = two*(x2 - x3);
      const Rational B2 = two*(y2 - y3);
      const Rational C2 = y3*y3 - y2*y2 + x3*x3 - x2*x2;

      // Compute the coordinates of the intersection point between the
      // two lines, given by (Nx / D, Ny / D), where:
      const Rational Nx = B1*C2 - B2*C1;
      const Rational Ny = A2*C1 - A1*C2;
      const Rational D = A1*B2 - A2*B1;

      // Make sure the three points are not collinear.
      const bool points_collinear = (CGAL::sign(D) == ZERO);

      if (points_collinear) {
        arc.reset_flags();            // inavlid arc
        return arc;
      }

      // The equation of the underlying circle is given by:
      Rational rat_coeffs[6];
      rat_coeffs[0] = D*D;
      rat_coeffs[1] = D*D;
      rat_coeffs[2] = 0;
      rat_coeffs[3] = -two*D*Nx;
      rat_coeffs[4] = -two*D*Ny;
      rat_coeffs[5] =
        Nx*Nx + Ny*Ny - ((D*x2 - Nx)*(D*x2 - Nx) + (D*y2 - Ny)*(D*y2 - Ny));

      // Determine the orientation: If the mid-point forms a left-turn with
      // the source and the target points, the orientation is positive (going
      // counterclockwise).
      // Otherwise, it is negative (going clockwise).
      auto orient_f = alg_kernel->orientation_2_object();
      Point_2 p_mid(nt_traits->convert(x2), nt_traits->convert(y2));

      auto orient = (orient_f(source, p_mid, target) == LEFT_TURN) ?
        COUNTERCLOCKWISE : CLOCKWISE;
      arc.set_orientation(orient);

      // Set the arc properties (no need to compute the orientation).
      m_traits.set(arc, rat_coeffs);
      return arc;
    }

    /*! Construct a conic arc from given five points, specified by the
     * points p1, p2, p3, p4 and p5.
     * \param p1 The source point of the given arc.
     * \param p2,p3,p4 Points lying on the conic arc, between p1 and p5.
     * \param p5 The target point of the given arc.
     * \pre No three points are collinear.
     */
    Curve_2 operator()(const Rat_point_2& p1, const Rat_point_2& p2,
                       const Rat_point_2& p3, const Rat_point_2& p4,
                       const Rat_point_2& p5) const {
      Curve_2 arc;

      // Make sure that no three points are collinear.
      auto orient_f = m_traits.m_rat_kernel->orientation_2_object();
      const bool point_collinear =
        (orient_f(p1, p2, p3) == COLLINEAR ||
         orient_f(p1, p2, p4) == COLLINEAR ||
         orient_f(p1, p2, p5) == COLLINEAR ||
         orient_f(p1, p3, p4) == COLLINEAR ||
         orient_f(p1, p3, p5) == COLLINEAR ||
         orient_f(p1, p4, p5) == COLLINEAR ||
         orient_f(p2, p3, p4) == COLLINEAR ||
         orient_f(p2, p3, p5) == COLLINEAR ||
         orient_f(p2, p4, p5) == COLLINEAR ||
         orient_f(p3, p4, p5) == COLLINEAR);

      if (point_collinear) {
        arc.reset_flags();            // inavlid arc
        return arc;
      }

      // Set the source and target.
      const Rational& x1 = p1.x();
      const Rational& y1 = p1.y();
      const Rational& x5 = p5.x();
      const Rational& y5 = p5.y();

      const auto nt_traits = m_traits.m_nt_traits;
      Point_2 source(nt_traits->convert(x1), nt_traits->convert(y1));
      Point_2 target(nt_traits->convert(x5), nt_traits->convert(y5));
      arc.set_endpoints(source, target);

      // Set a conic curve that passes through the five given point.
      typename Rat_kernel::Conic_2 temp_conic;
      temp_conic.set(p1, p2, p3, p4, p5);

      // Get the conic coefficients.
      Rational rat_coeffs[6];
      rat_coeffs[0] = temp_conic.r();
      rat_coeffs[1] = temp_conic.s();
      rat_coeffs[2] = temp_conic.t();
      rat_coeffs[3] = temp_conic.u();
      rat_coeffs[4] = temp_conic.v();
      rat_coeffs[5] = temp_conic.w();

      // Determine the orientation: If one of the midpoints forms a left-turn
      // with the source and the target points, the orientation is positive
      // (going counterclockwise).
      // Otherwise, it is negative (going clockwise).
      const Orientation turn = orient_f(p1, p2, p5);

      if (turn == LEFT_TURN) {
        arc.set_orientation(COUNTERCLOCKWISE);
        CGAL_precondition(orient_f(p1, p3, p5) == LEFT_TURN &&
                          orient_f(p1, p4, p5) == LEFT_TURN);
      }
      else {
        arc.set_orientation(CLOCKWISE);
        CGAL_precondition(orient_f(p1, p3, p5) != LEFT_TURN &&
                          orient_f(p1, p4, p5) != LEFT_TURN);
      }

      // Set the arc properties (no need to compute the orientation).
      m_traits.set(arc, rat_coeffs);

      // Make sure that all midpoints are strictly between the
      // source and the target.
      Point_2 mp2(nt_traits->convert(p2.x()), nt_traits->convert(p2.y()));
      Point_2 mp3(nt_traits->convert(p3.x()), nt_traits->convert(p3.y()));
      Point_2 mp4(nt_traits->convert(p4.x()), nt_traits->convert(p4.y()));

      if (! m_traits.is_strictly_between_endpoints(arc, mp2) ||
          ! m_traits.is_strictly_between_endpoints(arc, mp3) ||
          ! m_traits.is_strictly_between_endpoints(arc, mp4))
      {
        arc.reset_flags();            // inavlid arc
        return arc;
      }
      return arc;
    }

    /*! Construct a conic arc that lies on a conic given by its coefficients:
     *   `C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0`
     * The source and the target are specified by the intersection of the
     * conic with:
     *   `C_1: r_1*x^2 + s_1*y^2 + t_1*xy + u_1*x + v_1*y + w_1 = 0`
     *   `C_2: r_2*x^2 + s_2*y^2 + t_2*xy + u_2*x + v_2*y + w_2 = 0`
     * The user must also specify the source and the target with approximated
     * coordinates. The actual intersection points that best fits the source
     * (or the target) will be selected.
     */
    Curve_2 operator()(const Rational& r, const Rational& s, const Rational& t,
                       const Rational& u, const Rational& v, const Rational& w,
                       Orientation orient,
                       const Point_2& app_source,
                       const Rational& r_1, const Rational& s_1,
                       const Rational& t_1, const Rational& u_1,
                       const Rational& v_1, const Rational& w_1,
                       const Point_2& app_target,
                       const Rational& r_2, const Rational& s_2,
                       const Rational& t_2, const Rational& u_2,
                       const Rational& v_2, const Rational& w_2) const {
      Curve_2 arc;
      arc.set_orientation(orient);

      // Create the integer coefficients of the base conic.
      Rational rat_coeffs [6];
      rat_coeffs[0] = r;
      rat_coeffs[1] = s;
      rat_coeffs[2] = t;
      rat_coeffs[3] = u;
      rat_coeffs[4] = v;
      rat_coeffs[5] = w;

      const auto nt_traits = m_traits.m_nt_traits;
      Integer base_coeffs[6];
      nt_traits->convert_coefficients(rat_coeffs, rat_coeffs + 6, base_coeffs);

      int deg_base = (CGAL::sign(base_coeffs[0]) == ZERO &&
                      CGAL::sign(base_coeffs[1]) == ZERO &&
                      CGAL::sign(base_coeffs[2]) == ZERO) ? 1 : 2;

      // Compute the endpoints.
      Rational aux_rat_coeffs [6];
      Algebraic xs[4];
      Algebraic ys[4];
      Algebraic val;
      bool found;
      double curr_dist;
      double min_dist = -1;
      Integer aux_coeffs[6];
      for (int k = 1; k <= 2; ++k) {
        // Get the integer coefficients of the k-th auxiliary conic curve.
        aux_rat_coeffs[0] = (k == 1) ? r_1 : r_2;
        aux_rat_coeffs[1] = (k == 1) ? s_1 : s_2;
        aux_rat_coeffs[2] = (k == 1) ? t_1 : t_2;
        aux_rat_coeffs[3] = (k == 1) ? u_1 : u_2;
        aux_rat_coeffs[4] = (k == 1) ? v_1 : v_2;
        aux_rat_coeffs[5] = (k == 1) ? w_1 : w_2;

        nt_traits->convert_coefficients(aux_rat_coeffs, aux_rat_coeffs + 6,
                                        aux_coeffs);

        int deg_aux = ((CGAL::sign(aux_coeffs[0]) == ZERO) &&
                       (CGAL::sign(aux_coeffs[1]) == ZERO) &&
                       (CGAL::sign(aux_coeffs[2]) == ZERO)) ? 1 : 2;

        // Compute the x- and y-coordinates of intersection points of the base
        // conic and the k-th auxiliary conic.
        int n_xs = compute_resultant_roots(*nt_traits,
                                           base_coeffs[0], base_coeffs[1],
                                           base_coeffs[2],
                                           base_coeffs[3], base_coeffs[4],
                                           base_coeffs[5],
                                           deg_base,
                                           aux_coeffs[0], aux_coeffs[1],
                                           aux_coeffs[2],
                                           aux_coeffs[3], aux_coeffs[4],
                                           aux_coeffs[5],
                                           deg_aux,
                                           xs);

        int n_ys = compute_resultant_roots(*nt_traits,
                                           base_coeffs[1], base_coeffs[0],
                                           base_coeffs[2],
                                           base_coeffs[4], base_coeffs[3],
                                           base_coeffs[5],
                                           deg_base,
                                           aux_coeffs[1], aux_coeffs[0],
                                           aux_coeffs[2],
                                           aux_coeffs[4], aux_coeffs[3],
                                           aux_coeffs[5],
                                           deg_aux,
                                           ys);

        // Find the intersection point which is nearest the given approximation
        // and set it as the endpoint.
        found = false;
        for (int i = 0; i < n_xs; ++i) {
          for (int j = 0; j < n_ys; ++j) {
            // Check if the point (xs[i], ys[j]) lies on both conics.
            val = nt_traits->convert(base_coeffs[0]) * xs[i]*xs[i] +
              nt_traits->convert(base_coeffs[1]) * ys[j]*ys[j] +
              nt_traits->convert(base_coeffs[2]) * xs[i]*ys[j] +
              nt_traits->convert(base_coeffs[3]) * xs[i] +
              nt_traits->convert(base_coeffs[4]) * ys[j] +
              nt_traits->convert(base_coeffs[5]);

            if (CGAL::sign(val) != ZERO) continue;

            val = nt_traits->convert(aux_coeffs[0]) * xs[i]*xs[i] +
              nt_traits->convert(aux_coeffs[1]) * ys[j]*ys[j] +
              nt_traits->convert(aux_coeffs[2]) * xs[i]*ys[j] +
              nt_traits->convert(aux_coeffs[3]) * xs[i] +
              nt_traits->convert(aux_coeffs[4]) * ys[j] +
              nt_traits->convert(aux_coeffs[5]);

            if (CGAL::sign(val) == ZERO) {
              // Compute the distance of (xs[i], ys[j]) from the approximated
              // endpoint.
              double dx, dy;
              if (k == 1) {
                dx = CGAL::to_double (xs[i] - app_source.x());
                dy = CGAL::to_double (ys[j] - app_source.y());
              }
              else {
                dx = CGAL::to_double (xs[i] - app_target.x());
                dy = CGAL::to_double (ys[j] - app_target.y());
              }

              curr_dist = dx*dx + dy*dy;

              // Update the endpoint if (xs[i], ys[j]) is the nearest pair so
              // far.
              if (! found || curr_dist < min_dist) {
                if (k == 1) arc.set_source(Point_2(xs[i], ys[j]));
                else arc.set_target(Point_2(xs[i], ys[j]));
                min_dist = curr_dist;
                found = true;
              }
            }
          }
        }

        if (! found) {
          arc.reset_flags();          // inavlid arc
          return arc;
        }
      }

      // Make sure that the source and the target are not the same.
      auto cmp_xy = m_traits.m_alg_kernel->compare_xy_2_object();
      if (cmp_xy(arc.source(), arc.target()) == EQUAL) {
        arc.reset_flags();            // inavlid arc
        return arc;
      }

      // Set the arc properties (no need to compute the orientation).
      m_traits.set(arc, rat_coeffs);
      return arc;
    }

    /*! Return a segment connecting the two given endpoints.
     * \param source The source point.
     * \param target The target point.
     * \pre `source` and `target` must not be the same.
     * \return A segment connecting `source` and `target`.
     */
    Curve_2 operator()(const Point_2& source, const Point_2& target) const {
      const auto alg_kernel = m_traits.m_alg_kernel;
      CGAL_precondition_code(auto cmp_xy = alg_kernel->compare_xy_2_object());
      CGAL_precondition(cmp_xy(source, target) != EQUAL);

      Curve_2 cv;
      cv.set_coefficients(0, 0, 0, 0, 0, 0);
      cv.set_orientation(COLLINEAR);
      cv.set_flag(Curve_2::IS_VALID);
      cv.set_endpoints(source, target);
      cv.update_extra_data();
      return cv;
    }

    /*! Construct a conic arc from a given line segment.
     * \param seg The line segment with rational endpoints.
     */
    Curve_2 operator()(const Rat_segment_2& seg) const {
      Curve_2 cv;
      cv.set_orientation(COLLINEAR);

      // Set the source and target.
      const auto rat_kernel = m_traits.m_rat_kernel;
      Rat_point_2 source = rat_kernel->construct_vertex_2_object()(seg, 0);
      Rat_point_2 target = rat_kernel->construct_vertex_2_object()(seg, 1);
      const Rational& x1 = source.x();
      const Rational& y1 = source.y();
      const Rational& x2 = target.x();
      const Rational& y2 = target.y();

      const auto nt_traits = m_traits.m_nt_traits;
      cv.set_source(Point_2(nt_traits->convert(x1), nt_traits->convert(y1)));
      cv.set_target(Point_2(nt_traits->convert(x2), nt_traits->convert(y2)));

      // Make sure that the source and the taget are not the same.
      CGAL_precondition_code(auto cmp_xy = rat_kernel->compare_xy_2_object());
      CGAL_precondition(cmp_xy(source, target) != EQUAL);

      // The supporting conic is r=s=t=0, and u*x + v*y + w = 0 should hold
      // for both the source (x1,y1) and the target (x2, y2).
      const Rational zero(0);
      const Rational one(1);
      Rational rat_coeffs[6];

      rat_coeffs[0] = zero;
      rat_coeffs[1] = zero;
      rat_coeffs[2] = zero;

      if (rat_kernel->compare_x_2_object()(source, target) == EQUAL) {
        // The supporting conic is a vertical line, of the form x = CONST.
        rat_coeffs[3] = one;
        rat_coeffs[4] = zero;
        rat_coeffs[5] = -x1;
      }
      else {
        // The supporting line is A*x + B*y + C = 0, where:
        //
        //  A = y2 - y1,    B = x1 - x2,    C = x2*y1 - x1*y2
        //
        rat_coeffs[3] = y2 - y1;
        rat_coeffs[4] = x1 - x2;
        rat_coeffs[5] = x2*y1 - x1*y2;
      }

      // Set the arc properties (no need to compute the orientation).
      m_traits.set(cv, rat_coeffs);
      return cv;
    }

    /*! Construct a conic arc that is a full circle.
     * \param circ The circle with rational center and rational squared radius.
     */
    Curve_2 operator()(const Rat_circle_2& circ) const {
      Curve_2 cv;
      cv.set_orientation(CLOCKWISE);

      // Get the circle properties.
      const auto rat_kernel = m_traits.m_rat_kernel;
      Rat_point_2 center = rat_kernel->construct_center_2_object()(circ);
      Rational x0 = center.x();
      Rational y0 = center.y();
      Rational r_sqr = rat_kernel->compute_squared_radius_2_object()(circ);

      // Produce the correponding conic: if the circle center is (x0,y0)
      // and its squared radius is R^2, that its equation is:
      //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - R^2) = 0
      // Note that this equation describes a curve with a negative (clockwise)
      // orientation.
      const Rational zero(0);
      const Rational one(1);
      const Rational minus_two(-2);
      Rational rat_coeffs[6];
      rat_coeffs[0] = one;
      rat_coeffs[1] = one;
      rat_coeffs[2] = zero;
      rat_coeffs[3] = minus_two*x0;
      rat_coeffs[4] = minus_two*y0;
      rat_coeffs[5] = x0*x0 + y0*y0 - r_sqr;

      // Set the arc to be the full conic (no need to compute the orientation).
      m_traits.set_full(cv, rat_coeffs, false);
      return cv;
    }

    /*! Construct a conic arc that lies on a given circle:
     *   `C: (x - x0)^2 + (y - y0)^2 = R^2`
     * \param orient The orientation of the circle.
     * \param source The source point.
     * \param target The target point.
     * \pre The source and the target must be on the conic boundary and must
     *      not be the same.
     */
    Curve_2 operator()(const Rat_circle_2& circ, Orientation orient,
                       const Point_2& source, const Point_2& target) const {
      // Make sure that the source and the taget are not the same.
      CGAL_precondition_code(auto cmp_xy =
                               m_traits.m_alg_kernel->compare_xy_2_object());
      CGAL_precondition(cmp_xy(source, target) != EQUAL);
      CGAL_precondition(orient != COLLINEAR);

      Curve_2 cv;
      cv.set_endpoints(source, target);
      cv.set_orientation(orient);

      // Get the circle properties.
      const auto rat_kernel = m_traits.m_rat_kernel;
      Rat_point_2 center = rat_kernel->construct_center_2_object()(circ);
      Rational x0 = center.x();
      Rational y0 = center.y();
      Rational r_sqr = rat_kernel->compute_squared_radius_2_object()(circ);

      // Produce the correponding conic: if the circle center is (x0,y0)
      // and it squared radius is R^2, that its equation is:
      //   x^2 + y^2 - 2*x0*x - 2*y0*y + (x0^2 + y0^2 - R^2) = 0
      // Since this equation describes a curve with a negative (clockwise)
      // orientation, we multiply it by -1 if nece_Conic_arc_2 ssary to obtain a
      // positive (counterclockwise) orientation.
      const Rational zero(0);
      Rational rat_coeffs[6];

      if (cv.orientation() == COUNTERCLOCKWISE) {
        const Rational minus_one(-1);
        const Rational two(2);

        rat_coeffs[0] = minus_one;
        rat_coeffs[1] = minus_one;
        rat_coeffs[2] = zero;
        rat_coeffs[3] = two*x0;
        rat_coeffs[4] = two*y0;
        rat_coeffs[5] = r_sqr - x0*x0 - y0*y0;
      }
      else {
        const Rational one(1);
        const Rational minus_two(-2);

        rat_coeffs[0] = one;
        rat_coeffs[1] = one;
        rat_coeffs[2] = zero;
        rat_coeffs[3] = minus_two*x0;
        rat_coeffs[4] = minus_two*y0;
        rat_coeffs[5] = x0*x0 + y0*y0 - r_sqr;
      }

      // Set the arc properties (no need to compute the orientation).
      m_traits.set(cv, rat_coeffs);
      return cv;
    }
  };

  /*! Obtain a Construct_curve_2 functor object. */
  Construct_curve_2 construct_curve_2_object() const
  { return Construct_curve_2(*this); }
  //@}

  /// \name Functor definitions for the Boolean set-operation traits.
  //@{

  class Compare_endpoints_xy_2 {
  public:
    /*! Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param cv The curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv) const
    { return (cv.is_directed_right()) ? SMALLER : LARGER; }
  };

  /*! Obtain a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  class Construct_opposite_2 {
  public:
    /*! Construct an opposite x-monotone (with swapped source and target).
     * \param cv The curve.
     * \return The opposite curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& cv) const
    { return cv.flip(); }
  };

  /*! Obtain a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }

  class Trim_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Trim_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*!\brief
     * Returns a trimmed version of an cv
     *
     * \param xcv The arc
     * \param src the new first endpoint
     * \param tgt the new second endpoint
     * \return The trimmed arc
     *
     * \pre src != tgt
     * \pre both points must be interior and must lie on \c cv
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv,
                                  const Point_2& src, const Point_2& tgt) const
    {
      // make functor objects
      CGAL_precondition_code(Compare_y_at_x_2 compare_y_at_x_2 =
                             m_traits.compare_y_at_x_2_object());
      CGAL_precondition_code(Equal_2 equal_2 = m_traits.equal_2_object());
      Compare_x_2 compare_x_2 = m_traits.compare_x_2_object();
      // Check  whether source and taget are two distinct points and they lie
      // on the line.
      CGAL_precondition(compare_y_at_x_2(src, xcv) == EQUAL);
      CGAL_precondition(compare_y_at_x_2(tgt, xcv) == EQUAL);
      CGAL_precondition(! equal_2(src, tgt));

      //check if the orientation conforms to the src and tgt.
      return ((xcv.is_directed_right() && compare_x_2(src, tgt) == LARGER) ||
              (! xcv.is_directed_right() && compare_x_2(src, tgt) == SMALLER)) ?
        trim(xcv, tgt, src) : trim(xcv, src, tgt);
    }

  private:
    /*! Trim the arc given its new endpoints.
     * \param ps The new source point.
     * \param pt The new target point.
     * \return The new trimmed arc.
     * \pre Both ps and pt lies on the arc and must conform with the current
     *      direction of the arc.
    */
    X_monotone_curve_2 trim(const X_monotone_curve_2& xcv,
                            const Point_2& ps, const Point_2& pt) const {
      // Make sure that both ps and pt lie on the arc.
      CGAL_precondition(m_traits.contains_point(xcv, ps) &&
                        m_traits.contains_point(xcv, pt));

      X_monotone_curve_2 res_xcv = xcv; // make a copy of the current arc

      auto eq = m_traits.m_alg_kernel->equal_2_object();
      auto set_source = [&](const Point_2 ps)->void {
                          if (! eq(ps, xcv.source())) {
                            res_xcv.set_source(ps);
                            if (! ps.is_generating_conic(xcv.id()))
                              res_xcv.source().set_generating_conic(xcv.id());
                          }
                        };
      auto set_target = [&](const Point_2 pt)->void {
                          if (! eq(pt, xcv.target())) {
                            res_xcv.set_target(pt);
                            if (! pt.is_generating_conic(xcv.id()))
                              res_xcv.target().set_generating_conic(xcv.id());
                          }
                        };

      // Make sure that the endpoints conform with the direction of the arc.
      auto cmp_xy = m_traits.m_alg_kernel->compare_xy_2_object();
      auto res = cmp_xy(ps, pt);
      CGAL_assertion(res != EQUAL);
      if ((xcv.test_flag(X_monotone_curve_2::IS_DIRECTED_RIGHT) &&
           (res == LARGER)) ||
          (! xcv.test_flag(X_monotone_curve_2::IS_DIRECTED_RIGHT) &&
           (res == SMALLER))) {
        set_source(pt);
        set_target(ps);
      }
      else {
        set_source(ps);
        set_target(pt);
      }

      return res_xcv;
    }
  };

  /*! Obtain a Trim_2 functor object. */
  Trim_2 trim_2_object() const { return Trim_2(*this); }
  //@}

  /// \name Extra functor definitions.
  //@{

  class Construct_bbox_2 {
  protected:
    using Traits = Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits The traits.
     */
    Construct_bbox_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Obtain a bounding box for the conic arc.
     * \return The bounding box.
     */
    Bbox_2 operator()(const X_monotone_curve_2& xcv) const { return bbox(xcv); }
    Bbox_2 operator()(const Curve_2& cv) const { return bbox(cv); }

  private:
    Bbox_2 bbox(const X_monotone_curve_2& xcv) const {
      CGAL_precondition(xcv.is_valid());
      double x_min(0), y_min(0), x_max(0), y_max(0);

      if (xcv.is_full_conic()) {
        // In case of a full conic (an ellipse or a circle), compute the
        // horizontal and vertical tangency points and use them to bound the arc.
        Alg_point_2 tan_ps[2];
        CGAL_assertion_code(size_t n_tan_ps);

        CGAL_assertion_code
          (n_tan_ps = m_traits.vertical_tangency_points(xcv, tan_ps));
        CGAL_assertion(n_tan_ps == 2);

        if (CGAL::to_double(tan_ps[0].x()) < CGAL::to_double(tan_ps[1].x())) {
          x_min = CGAL::to_double(tan_ps[0].x());
          x_max = CGAL::to_double(tan_ps[1].x());
        }
        else {
          x_min = CGAL::to_double(tan_ps[1].x());
          x_max = CGAL::to_double(tan_ps[0].x());
        }

        CGAL_assertion_code
          (n_tan_ps = m_traits.horizontal_tangency_points(xcv, tan_ps));
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
          CGAL::to_double(xcv.source().x()) < CGAL::to_double(xcv.target().x());
        x_min = (source_left) ?
          CGAL::to_double(xcv.source().x()) :
          CGAL::to_double(xcv.target().x());
        x_max = source_left ?
          CGAL::to_double(xcv.target().x()) :
          CGAL::to_double(xcv.source().x());

        bool source_down =
          CGAL::to_double(xcv.source().y()) < CGAL::to_double(xcv.target().y());
        y_min = source_down ?
          CGAL::to_double(xcv.source().y()) :
          CGAL::to_double(xcv.target().y());
        y_max = source_down ?
          CGAL::to_double(xcv.target().y()) :
          CGAL::to_double(xcv.source().y());

        // Go over the vertical tangency points and try to update the x-points.
        Alg_point_2 tan_ps[2];
        auto n_tan_ps = m_traits.vertical_tangency_points(xcv, tan_ps);
        for (decltype(n_tan_ps) i = 0; i < n_tan_ps; ++i) {
          if (CGAL::to_double(tan_ps[i].x()) < x_min)
            x_min = CGAL::to_double(tan_ps[i].x());
          if (CGAL::to_double(tan_ps[i].x()) > x_max)
            x_max = CGAL::to_double(tan_ps[i].x());
        }

        // Go over the horizontal tangency points and try to update the y-points.
        n_tan_ps = m_traits.horizontal_tangency_points(xcv, tan_ps);
        for (decltype(n_tan_ps) i = 0; i < n_tan_ps; ++i) {
          if (CGAL::to_double(tan_ps[i].y()) < y_min)
            y_min = CGAL::to_double(tan_ps[i].y());
          if (CGAL::to_double(tan_ps[i].y()) > y_max)
            y_max = CGAL::to_double(tan_ps[i].y());
        }
      }

      // Return the resulting bounding box.
      return Bbox_2(x_min, y_min, x_max, y_max);
    }
  };

  /*! Obtain a Bbox_2 functor object. */
  Construct_bbox_2 construct_bbox_2_object() const
  { return Construct_bbox_2(*this); }
  //@}

  /*! Set the properties of a conic arc (for the usage of the constructors).
   * \param rat_coeffs A vector of size 6, storing the rational coefficients
   *                   of \f$x^2\f$, \f$y^2\f$, \f$x \cdot y\f$, \f$x\f$, \f$y\f$
   *                   and the free coefficient resp.
   */
  void set(Curve_2& cv, const Rational* rat_coeffs) const {
    cv.set_flag(Curve_2::IS_VALID);

    // Convert the coefficients vector to an equivalent vector of integer
    // coefficients.
    Integer int_coeffs[6];
    m_nt_traits->convert_coefficients(rat_coeffs, rat_coeffs + 6, int_coeffs);

    // Check the orientation of conic curve, and negate the conic coefficients
    // if its given orientation.
    typename Rat_kernel::Conic_2 temp_conic(rat_coeffs[0], rat_coeffs[1],
                                            rat_coeffs[2], rat_coeffs[3],
                                            rat_coeffs[4], rat_coeffs[5]);


    Integer r, s, t, u, v, w;
    if (cv.orientation() == temp_conic.orientation()) {
      r = int_coeffs[0];
      s = int_coeffs[1];
      t = int_coeffs[2];
      u = int_coeffs[3];
      v = int_coeffs[4];
      w = int_coeffs[5];
    }
    else {
      r = -int_coeffs[0];
      s = -int_coeffs[1];
      t = -int_coeffs[2];
      u = -int_coeffs[3];
      v = -int_coeffs[4];
      w = -int_coeffs[5];
    }
    cv.set_coefficients(r, s, t, u, v, w);

    const auto& source = cv.source();
    const auto& target = cv.target();
    // Make sure both endpoint lie on the supporting conic.
    if (! is_on_supporting_conic(cv, source) ||
        ! is_on_supporting_conic(cv, target))
    {
      cv.reset_flags();            // inavlid arc
      return;
    }

    // Check whether we have a degree 2 curve.
    if ((CGAL::sign(r) != ZERO) || (CGAL::sign(s) != ZERO) ||
        (CGAL::sign(t) != ZERO))
    {
      if (cv.orientation() == COLLINEAR) {
        // Make sure the midpoint is on the line pair (thus making sure that
        // the two points are not taken from different lines).
        auto ctr_mid_point = m_alg_kernel->construct_midpoint_2_object();
        Point_2 p_mid = ctr_mid_point(source, target);
        // if (! is_on_supporting_conic(arc, p_mid))
        if (CGAL::sign((m_nt_traits->convert(r)*p_mid.x() +
                        m_nt_traits->convert(t)*p_mid.y() +
                        m_nt_traits->convert(u)) * p_mid.x() +
                       (m_nt_traits->convert(s)*p_mid.y() +
                        m_nt_traits->convert(v)) * p_mid.y() +
                       m_nt_traits->convert(w)) != ZERO)
        {
          cv.reset_flags();    // inavlid arc
          return;
        }

        // We have a segment of a line pair with rational coefficients.
        // Compose the equation of the underlying line
        // (with algebraic coefficients).
        cv.update_extra_data();
      }
      else {
        // The sign of (4rs - t^2) detetmines the conic type:
        // - if it is possitive, the conic is an ellipse,
        // - if it is negative, the conic is a hyperbola,
        // - if it is zero, the conic is a parabola.
        CGAL::Sign sign_conic = CGAL::sign(4*r*s - t*t);

        // Build the extra hyperbolic data if necessary
        if (sign_conic == NEGATIVE) build_hyperbolic_arc_data(cv);

        if (sign_conic != POSITIVE) {
          // In case of a non-degenerate parabola or a hyperbola, make sure
          // the arc is not infinite.
          Point_2 p_mid =
            m_alg_kernel->construct_midpoint_2_object()(source, target);
          Alg_point_2 ps[2];
          bool finite_at_x = (points_at_x(cv, p_mid, ps) > 0);
          bool finite_at_y = (points_at_y(cv, p_mid, ps) > 0);
          if (! finite_at_x && ! finite_at_y) {
            cv.reset_flags();              // inavlid arc
            return;
          }
        }
      }
    }


    cv.set_flag(Curve_2::IS_VALID);            // arc is valid
    cv.reset_flag(Curve_2::IS_FULL_CONIC);     // not a full conic
  }

  /*! Set the properties of a conic arc that is really a full curve
   * (that is, an ellipse).
   * \param rat_coeffs A vector of size 6, storing the rational coefficients
   *                   of \f$x^2\f$, \f$y^2\f$, \f$x \cdot y\f$, \f$x\f$, \f$y\f$
   *                   and the free coefficient resp.
   * \param comp_orient Should we compute the orientation of the given curve.
   */
  void set_full(Curve_2& cv, const Rational* rat_coeffs,
                const bool& comp_orient) const {
    // Convert the coefficients vector to an equivalent vector of integer
    // coefficients.
    Integer int_coeffs[6];
    m_nt_traits->convert_coefficients(rat_coeffs, rat_coeffs + 6, int_coeffs);

    // Check the orientation of conic curve, and negate the conic coefficients
    // if its given orientation.
    typename Rat_kernel::Conic_2 temp_conic(rat_coeffs[0], rat_coeffs[1],
                                            rat_coeffs[2], rat_coeffs[3],
                                            rat_coeffs[4], rat_coeffs[5]);
    const Orientation temp_orient = temp_conic.orientation();

    if (comp_orient) cv.set_orientation(temp_orient);

    Integer r, s, t, u, v, w;
    if (cv.orientation() == temp_orient) {
      r = int_coeffs[0];
      s = int_coeffs[1];
      t = int_coeffs[2];
      u = int_coeffs[3];
      v = int_coeffs[4];
      w = int_coeffs[5];
    }
    else {
      r = -int_coeffs[0];
      s = -int_coeffs[1];
      t = -int_coeffs[2];
      u = -int_coeffs[3];
      v = -int_coeffs[4];
      w = -int_coeffs[5];
    }
    cv.set_coefficients(r, s, t, u, v, w);

    // Make sure the conic is a non-degenerate ellipse:
    // The coefficients should satisfy (4rs - t^2) > 0.
    const bool is_ellipse = (CGAL::sign(4*r*s - t*t) == POSITIVE);
    CGAL_assertion(is_ellipse);

    // We do not have to store any extra data with the arc.

    // Mark that this arc is a full conic curve.
    if (is_ellipse) {
      cv.set_flag(Curve_2::IS_VALID);
      cv.set_flag(Curve_2::IS_FULL_CONIC);
    }
    else cv.reset_flags();            // inavlid arc
  }

  /*! Check whether the given point lies on the supporting conic of the arc.
   * \param p The query point.
   * \return true if p lies on the supporting conic; (false) otherwise.
   */
  bool is_on_supporting_conic(Curve_2& cv, const Point_2& p) const {
    // Check whether p satisfies the conic equation.
    const Algebraic r = m_nt_traits->convert(cv.r());
    const Algebraic t = m_nt_traits->convert(cv.t());
    const Algebraic u = m_nt_traits->convert(cv.u());
    const Algebraic s = m_nt_traits->convert(cv.s());
    const Algebraic v = m_nt_traits->convert(cv.v());
    const Algebraic w = m_nt_traits->convert(cv.w());

    // The point must satisfy: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0.
    const Algebraic val =
      (r*p.x() + t*p.y() + u)*p.x() + (s*p.y() + v)* p.y() + w;
    return (CGAL::sign(val) == ZERO);
  }

  /*! Check whether the given point is between the source and the target.
   * The point is assumed to be on the conic's boundary.
   * \param p The query point.
   * \return `true` if the point is between the two endpoints;
   *         `false` if it is not.
   */
  bool is_between_endpoints(const Curve_2& cv, const Point_2& p) const {
    CGAL_precondition(! cv.is_full_conic());

    // Check if p is one of the endpoints.
    auto eq = m_alg_kernel->equal_2_object();
    if (eq(p, cv.source()) || eq(p, cv.target())) return true;
    else return is_strictly_between_endpoints(cv, p);
  }

  /*! Check whether the given point is strictly between the source and the
   * target (but not any of them).
   * The point is assumed to be on the conic's boundary.
   * \param p The query point.
   * \return `true` if the point is strictly between the two endpoints;
   *         `false` if it is not.
   */
  bool is_strictly_between_endpoints(const Curve_2& cv, const Point_2& p) const
  {
    // In case this is a full conic, any point on its boundary is between
    // its end points.
    if (cv.is_full_conic()) return true;

    // Check if we have extra data available.
    const auto* extra_data = cv.extra_data();
    if (extra_data != nullptr) {
      if (extra_data->side != ZERO) {
        // In case of a hyperbolic arc, make sure the point is located on the
        // same branch as the arc.
        if (cv.sign_of_extra_data(p.x(), p.y()) != extra_data->side)
          return false;
      }
      else {
        // In case we have a segment of a line pair, make sure that p really
        // satisfies the equation of the line.
        if (cv.sign_of_extra_data(p.x(), p.y()) != ZERO) return false;
      }
    }

    // Act according to the conic degree.
    auto orient_f = m_alg_kernel->orientation_2_object();
    const auto& source = cv.source();
    const auto& target = cv.target();
    if (cv.orientation() == COLLINEAR) {
      Comparison_result res1;
      Comparison_result res2;

      if (m_alg_kernel->compare_x_2_object()(source, target) == EQUAL) {
        // In case of a vertical segment - just check whether the y coordinate
        // of p is between those of the source's and of the target's.
        auto cmp_y = m_alg_kernel->compare_y_2_object();
        res1 = cmp_y(p, source);
        res2 = cmp_y(p, target);
      }
      else {
        // Otherwise, since the segment is x-monotone, just check whether the
        // x coordinate of p is between those of the source's and of the
        // target's.
        auto cmp_x = m_alg_kernel->compare_x_2_object();
        res1 = cmp_x(p, source);
        res2 = cmp_x(p, target);
      }

      // If p is not in the (open) x-range (or y-range) of the segment, it
      // cannot be contained in the segment.
      if ((res1 == EQUAL) || (res2 == EQUAL) || (res1 == res2)) return false;

      // Perform an orientation test: This is crucial for segment of line
      // pairs, as we want to make sure that p lies on the same line as the
      // source and the target.
      return (orient_f(source, p, target) == COLLINEAR);
    }

    // In case of a conic of degree 2, make a decision based on the conic's
    // orientation and whether (source,p,target) is a right or a left turn.
    if (cv.orientation() == COUNTERCLOCKWISE)
      return (orient_f(source, p, target) == LEFT_TURN);
    else
      return (orient_f(source, p, target) == RIGHT_TURN);
  }

  /*! Build the data for hyperbolic arc, contaning the characterization of the
   * hyperbolic branch the arc is placed on.
   */
  void build_hyperbolic_arc_data(Curve_2& cv) const {
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
    const int or_fact = (cv.orientation() == CLOCKWISE) ? -1 : 1;
    const Algebraic r = m_nt_traits->convert(Integer(or_fact * cv.r()));
    const Algebraic s = m_nt_traits->convert(Integer(or_fact * cv.s()));
    const Algebraic t = m_nt_traits->convert(Integer(or_fact * cv.t()));
    const Algebraic cos_2phi = (r - s) / m_nt_traits->sqrt((r-s)*(r-s) + t*t);
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
      sin_phi = m_nt_traits->sqrt((one + cos_2phi) / two);
      cos_phi = m_nt_traits->sqrt((one - cos_2phi) / two);
    }
    else {
      // sin(2*phi) < 0 so PI/2 < phi < PI.
      sin_phi = m_nt_traits->sqrt((one + cos_2phi) / two);
      cos_phi = - m_nt_traits->sqrt((one - cos_2phi) / two);
    }

    // Calculate the center (x0, y0) of the conic, given by the formulae:
    //
    //        t*v - 2*s*u                t*u - 2*r*v
    //  x0 = -------------   ,     y0 = -------------
    //        4*r*s - t^2                4*r*s - t^2
    //
    // The denominator (4*r*s - t^2) must be negative for hyperbolas.
    const Algebraic u = m_nt_traits->convert(Integer(or_fact * cv.u()));
    const Algebraic v = m_nt_traits->convert(Integer(or_fact * cv.v()));
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
    // We use it to make sure that the two endpoints are located on the same
    // branch of the hyperbola.
    auto a = cos_phi;
    auto b = sin_phi;
    auto c = - (cos_phi*x0 + sin_phi*y0);
    const auto& source = cv.source();
    auto val = a * source.x() + b * source.y() + c;
    auto side = CGAL::sign(val);
    CGAL_assertion(side != ZERO);
    cv.set_extra_data(a, b, c, side);
    CGAL_assertion_code(const auto& target = cv.target());
    CGAL_assertion(side == cv.sign_of_extra_data(target.x(), target.y()));
  }

  /*! Find the x coordinates of the underlying conic at a given y coordinate.
   * \param y The y coordinate.
   * \param xs The output x coordinates.
   * \pre The vector xs must be allocated at the size of 2.
   * \return The number of x coordinates computed (either 0, 1 or 2).
   */
  int conic_get_x_coordinates(const Curve_2& cv,
                              const Algebraic& y, Algebraic* xs) const {
    // Solve the quadratic equation for a given y and find the x values:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    Algebraic A = m_nt_traits->convert(cv.r());
    Algebraic B = m_nt_traits->convert(cv.t())*y + m_nt_traits->convert(cv.u());
    Algebraic C =
      (m_nt_traits->convert(cv.s())*y + m_nt_traits->convert(cv.v()))*y +
      m_nt_traits->convert(cv.w());

    return solve_quadratic_equation(A, B, C, xs[0], xs[1]);
  }

  /*! Find the y coordinates of the underlying conic at a given x coordinate.
   * \param x The x coordinate.
   * \param ys The output y coordinates.
   * \pre The vector ys must be allocated at the size of 2.
   * \return The number of y coordinates computed (either 0, 1 or 2).
   */
  int conic_get_y_coordinates(const Curve_2& cv,
                              const Algebraic& x, Algebraic* ys) const {
    // Solve the quadratic equation for a given x and find the y values:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    Algebraic A = m_nt_traits->convert(cv.s());
    Algebraic B = m_nt_traits->convert(cv.t())*x + m_nt_traits->convert(cv.v());
    Algebraic C =
      (m_nt_traits->convert(cv.r())*x + m_nt_traits->convert(cv.u()))*x +
      m_nt_traits->convert(cv.w());

    return solve_quadratic_equation(A, B, C, ys[0], ys[1]);
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
    Algebraic sqrt_disc = m_nt_traits->sqrt(disc);

    x_minus = -(B + sqrt_disc) / _2A;
    x_plus = (sqrt_disc - B) / _2A;
    return 2;
  }

  /*! Compute a point on an arc with the same \f$x\f$-coordiante as the given
   * point.
   * \param p The given point.
   * \pre The arc is not vertical and `p` is in the \f$x\f$-range of the arc.
   * \return A point on the arc with the same \f$x\f$-coordiante as `p`.
   */
  Point_2 point_at_x(const X_monotone_curve_2& xcv, const Point_2& p) const {
    // Make sure that p is in the x-range of the arc.
    CGAL_precondition(! xcv.is_vertical());

    CGAL_precondition_code(auto cmp_x = m_alg_kernel->compare_x_2_object());
    CGAL_precondition((cmp_x(p, xcv.left()) != SMALLER) &&
                      (cmp_x(p, xcv.right()) != LARGER));

    if (xcv.is_special_segment()) {
      // In case of a special segment, the equation of the supported line
      // (a*x + b*y + c) = 0 is stored with the extra data field, and we
      // simply have:
      const auto& extra_data = xcv.extra_data();
      Algebraic y = -(extra_data->a*p.x() + extra_data->c) / extra_data->b;

      // Return the computed point.
      return Point_2(p.x(), y);
    }

    // Compute the y-coordinate according to the degree of the supporting
    // conic curve.
    Algebraic y;

    if (xcv.degree_mask() == X_monotone_curve_2::degree_1_mask()) {
      // In case of a linear curve, the y-coordinate is a simple linear
      // expression of x(p) (note that v is not 0 as the arc is not vertical):
      //   y = -(u*x(p) + w) / v
      y = -(xcv.alg_u()*p.x() + xcv.alg_w()) / xcv.alg_v();
    }
    else if (xcv.orientation() == COLLINEAR) {
      const auto& extra_data = xcv.extra_data();
      CGAL_assertion(extra_data != nullptr);

      // In this case the equation of the supporting line is given by the
      // extra data structure.
      y = -(extra_data->a * p.x() + extra_data->c) / extra_data->b;
    }
    else {
      CGAL_assertion(xcv.degree_mask() == X_monotone_curve_2::degree_2_mask());

      // In this case the y-coordinate is one of solutions to the quadratic
      // equation:
      //  s*y^2 + (t*x(p) + v)*y + (r*x(p)^2 + u*x(p) + w) = 0
      Algebraic A = xcv.alg_s();
      Algebraic B = xcv.alg_t()*p.x() + xcv.alg_v();
      Algebraic C = (xcv.alg_r()*p.x() + xcv.alg_u())*p.x() + xcv.alg_w();

      if (CGAL::sign(xcv.s()) == ZERO) {
        // In this case A is 0 and we have a linear equation.
        CGAL_assertion(CGAL::sign(B) != ZERO);

        y = -C / B;
      }
      else {
        // Solve the quadratic equation.
        Algebraic disc = B*B - 4*A*C;

        CGAL_assertion(CGAL::sign(disc) != NEGATIVE);

        // We take either the root involving -sqrt(disc) or +sqrt(disc)
        // based on the information flags.
        y = (xcv.test_flag(X_monotone_curve_2::PLUS_SQRT_DISC_ROOT)) ?
          (m_nt_traits->sqrt(disc) - B) / (2*A) :
          -(B + m_nt_traits->sqrt(disc)) / (2*A);
      }
    }

    // Return the computed point.
    return Point_2(p.x(), y);
  }

  /*! Find all points on the arc with a given \f$x\f$-coordinate.
   * \param p A placeholder for the \f$x\f$-coordinate.
   * \param ps The point on the arc at `x(p)`.
   * \pre The vector `ps` should be allocated at the size of 2.
   * \return The number of points found.
   */
  int points_at_x(const Curve_2& cv, const Point_2& p, Alg_point_2* ps) const {
    // Get the y coordinates of the points on the conic.
    Algebraic ys[2];
    int n = conic_get_y_coordinates(cv, p.x(), ys);

    // Find all the points that are contained in the arc.
    int m = 0;

    for (int i = 0; i < n; ++i) {
      ps[m] = Point_2(p.x(), ys[i]);
      if (cv.is_full_conic() || is_between_endpoints(cv, ps[m])) ++m;
    }

    // Return the number of points on the arc.
    CGAL_assertion(m <= 2);
    return m;
  }

  /*! Find all points on the arc with a given \f$y\f$-coordinate.
   * \param p A placeholder for the \f$y\f$-coordinate.
   * \param ps The point on the arc at `x(p)`.
   * \pre The vector `ps` should be allocated at the size of 2.
   * \return The number of points found.
   */
  int points_at_y(const Curve_2& cv, const Point_2& p, Alg_point_2* ps) const {
    // Get the y coordinates of the points on the conic.
    Algebraic xs[2];
    int n = conic_get_x_coordinates(cv, p.y(), xs);

    // Find all the points that are contained in the arc.
    int m = 0;

    for (int i = 0; i < n; ++i) {
      ps[m] = Point_2(xs[i], p.y());
      if (cv.is_full_conic() || is_between_endpoints(cv, ps[m])) ++m;
    }

    // Return the number of points on the arc.
    CGAL_assertion(m <= 2);
    return m;
  }

  /*! Set the properties of the x-monotone conic arc (for the usage of the
   * constructors).
   */
  void set_x_monotone(X_monotone_curve_2& xcv) const {
    // Convert the coefficients of the supporting conic to algebraic numbers.
    auto alg_r = m_nt_traits->convert(xcv.r());
    auto alg_s = m_nt_traits->convert(xcv.s());
    auto alg_t = m_nt_traits->convert(xcv.t());
    auto alg_u = m_nt_traits->convert(xcv.u());
    auto alg_v = m_nt_traits->convert(xcv.v());
    auto alg_w = m_nt_traits->convert(xcv.w());
    xcv.set_alg_coefficients(alg_r, alg_s, alg_t, alg_u, alg_v, alg_w);

    // Set the generating conic ID for the source and target points.
    xcv.set_generating_conic(xcv.id());

    // Update the m_info bits.
    xcv.set_flag(Curve_2::IS_VALID);
    xcv.reset_flag(Curve_2::IS_FULL_CONIC);

    // Check (i) whether the arc is a vertical segment, and (ii) whether it is
    // directed right.
    auto cmp_x = m_alg_kernel->compare_x_2_object();
    Comparison_result resx = cmp_x(xcv.source(), xcv.target());
    if (resx == SMALLER) xcv.set_flag(X_monotone_curve_2::IS_DIRECTED_RIGHT);
    else if (resx == EQUAL) {
      xcv.set_flag(X_monotone_curve_2::IS_VERTICAL_SEGMENT);

      auto cmp_y = m_alg_kernel->compare_y_2_object();
      Comparison_result resy =  cmp_y(xcv.source(), xcv.target());
      CGAL_assertion(resy != EQUAL);
      if (resy == SMALLER) xcv.set_flag(X_monotone_curve_2::IS_DIRECTED_RIGHT);
    }

    if (xcv.orientation() == COLLINEAR) {
      // Compute the degree of the underlying conic.
      if ((CGAL::sign(xcv.r()) != ZERO) ||
          (CGAL::sign(xcv.s()) != ZERO) ||
          (CGAL::sign(xcv.t()) != ZERO))
      {
        xcv.set_flag(X_monotone_curve_2::DEGREE_2);
        xcv.set_flag(X_monotone_curve_2::IS_SPECIAL_SEGMENT);
      }
      else {
        xcv.set_flag(X_monotone_curve_2::DEGREE_1);

        // Check whether this is a special segment
        if ((CGAL::sign(xcv.u()) == ZERO) && (CGAL::sign(xcv.v()) == ZERO)) {
          xcv.set_flag(X_monotone_curve_2::IS_SPECIAL_SEGMENT);
          xcv.update_extra_data();
        }
      }
      return;
    }

    xcv.set_flag(X_monotone_curve_2::DEGREE_2);

    // Compute a midpoint between the source and the target and get the y-value
    // of the arc at its x-coordiante.
    Point_2 p_mid =
      m_alg_kernel->construct_midpoint_2_object()(xcv.source(), xcv.target());
    Algebraic ys[2];
    CGAL_assertion_code(int n_ys = )
      conic_get_y_coordinates(xcv, p_mid.x(), ys);

    CGAL_assertion(n_ys != 0);

    // Check which solution lies on the x-monotone arc.
    Point_2 p_arc_mid(p_mid.x(), ys[0]);

    if (is_strictly_between_endpoints(xcv, p_arc_mid)) {
      // Mark that we should use the -sqrt(disc) root for points on this
      // x-monotone arc.
      xcv.reset_flag(X_monotone_curve_2::PLUS_SQRT_DISC_ROOT);
    }
    else {
      CGAL_assertion(n_ys == 2);
      p_arc_mid = Point_2(p_mid.x(), ys[1]);
      CGAL_assertion(is_strictly_between_endpoints(xcv, p_arc_mid));

      // Mark that we should use the +sqrt(disc) root for points on this
      // x-monotone arc.
      xcv.set_flag(X_monotone_curve_2::PLUS_SQRT_DISC_ROOT);
    }

    // Check whether the conic is facing up or facing down:
    // Check whether the arc (which is x-monotone of degree 2) lies above or
    // below the segement that connects its two end-points (x1,y1) and (x2,y2).
    // To do that, we find the y coordinate of a point on the arc whose x
    // coordinate is (x1+x2)/2 and compare it to (y1+y2)/2.
    auto cmp_y = m_alg_kernel->compare_y_2_object();
    Comparison_result res = cmp_y(p_arc_mid, p_mid);

    // If the arc is above the connecting segment, so it is facing upwards.
    if (res == LARGER) xcv.set_flag(X_monotone_curve_2::FACING_UP);
    // If the arc is below the connecting segment, so it is facing downwards.
    else if (res == SMALLER) xcv.set_flag(X_monotone_curve_2::FACING_DOWN);
  }

  /*! Check whether the two arcs have the same supporting conic.
   * \param xcv1 The first comparedb arc.
   * \param xcv2 The secind compared arc.
   * \return `true` if the two supporting conics are the same.
   */
  bool has_same_supporting_conic(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2) const {
    // Check if the two arcs originate from the same conic:
    if ((xcv1.id() == xcv2.id()) &&
        xcv1.id().is_valid() && xcv2.id().is_valid())
      return true;

    // In case both arcs are collinear, check if they have the same
    // supporting lines.
    if ((xcv1.orientation() == COLLINEAR) &&
        (xcv2.orientation() == COLLINEAR)) {
      // Construct the two supporting lines and compare them.
      auto ctr_line = m_alg_kernel->construct_line_2_object();
      typename Alg_kernel::Line_2 l1 = ctr_line(xcv1.source(), xcv1.target());
      typename Alg_kernel::Line_2 l2 = ctr_line(xcv2.source(), xcv2.target());
      auto equal = m_alg_kernel->equal_2_object();
      if (equal(l1, l2)) return true;

      l2 = ctr_line(xcv2.target(), xcv2.source());
      return equal(l1, l2);
    }
    else if ((xcv1.orientation() == COLLINEAR) ||
             (xcv2.orientation() == COLLINEAR)) {
      // Only one arc is collinear; therefore so the supporting curves cannot
      // be the same:
      return false;
    }

    // Check whether the coefficients of the two supporting conics are equal
    // up to a constant factor.
    Integer factor1 = 1;
    Integer factor2 = 1;

    if (CGAL::sign(xcv1.r()) != ZERO) factor1 = xcv1.r();
    else if (CGAL::sign(xcv1.s()) != ZERO) factor1 = xcv1.s();
    else if (CGAL::sign(xcv1.t()) != ZERO) factor1 = xcv1.t();
    else if (CGAL::sign(xcv1.u()) != ZERO) factor1 = xcv1.u();
    else if (CGAL::sign(xcv1.v()) != ZERO) factor1 = xcv1.v();
    else if (CGAL::sign(xcv1.w()) != ZERO) factor1 = xcv1.w();

    if (CGAL::sign(xcv2.r()) != ZERO) factor2 = xcv2.r();
    else if (CGAL::sign(xcv2.s()) != ZERO) factor2 = xcv2.s();
    else if (CGAL::sign(xcv2.t()) != ZERO) factor2 = xcv2.t();
    else if (CGAL::sign(xcv2.u()) != ZERO) factor2 = xcv2.u();
    else if (CGAL::sign(xcv2.v()) != ZERO) factor2 = xcv2.v();
    else if (CGAL::sign(xcv2.w()) != ZERO) factor2 = xcv2.w();

    return (CGAL::compare(xcv1.r() * factor2, xcv2.r() * factor1) == EQUAL &&
            CGAL::compare(xcv1.s() * factor2, xcv2.s() * factor1) == EQUAL &&
            CGAL::compare(xcv1.t() * factor2, xcv2.t() * factor1) == EQUAL &&
            CGAL::compare(xcv1.u() * factor2, xcv2.u() * factor1) == EQUAL &&
            CGAL::compare(xcv1.v() * factor2, xcv2.v() * factor1) == EQUAL &&
            CGAL::compare(xcv1.w() * factor2, xcv2.w() * factor1) == EQUAL);
  }

  /*! Check whether the given point lies on the arc.
   * \param p The qury point.
   * \param (true) if p lies on the arc; (false) otherwise.
   */
  bool contains_point(const X_monotone_curve_2& xcv, const Point_2& p) const {
    // First check if p lies on the supporting conic. We first check whether
    // it is one of p's generating conic curves.
    bool p_on_conic(false);
    if (p.is_generating_conic(xcv.id())) p_on_conic = true;
    else {
      // Check whether p satisfies the supporting conic equation.
      p_on_conic = xcv.is_on_supporting_conic(p.x(), p.y());
      if (p_on_conic) {
        // As p lies on the supporting conic of our arc, add its ID to
        // the list of generating conics for p.
        Point_2& p_non_const = const_cast<Point_2&>(p);
        p_non_const.set_generating_conic(xcv.id());
      }
    }

    if (! p_on_conic) return false;

    // Check if p is between the endpoints of the arc.
    return is_between_endpoints(xcv, p);
  }

  /*! Find the vertical tangency points of the undelying conic.
   * \param ps The output points of vertical tangency.
   *           This area must be allocated at the size of 2.
   * \return The number of vertical tangency points.
   */
  int conic_vertical_tangency_points(const Curve_2& cv, Alg_point_2* ps) const {
    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    if (CGAL::sign(cv.s()) == ZERO) return 0;

    // We are interested in the x coordinates where the quadratic equation:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    // has a single solution (obviously if s = 0, there are no such points).
    // We therefore demand that the discriminant of this equation is zero:
    //  (t*x + v)^2 - 4*s*(r*x^2 + u*x + w) = 0
    const Integer two(2);
    const Integer four(4);
    Algebraic xs[2];

    auto r = cv.r();
    auto s = cv.s();
    auto t = cv.t();
    auto u = cv.u();
    auto v = cv.v();
    auto w = cv.w();
    Algebraic* xs_end = m_nt_traits->solve_quadratic_equation(Integer(t*t - four*r*s),
                                                              Integer(two*t*v - four*s*u),
                                                              Integer(v*v - four*s*w),
                                                              xs);
    auto n_xs = static_cast<int>(xs_end - xs);

    // Find the y-coordinates of the vertical tangency points.
    Algebraic ys[2];
    Algebraic* ys_end;
    int n_ys;

    if (CGAL::sign(cv.t()) == ZERO) {
      // The two vertical tangency points have the same y coordinate:
      ys[0] = m_nt_traits->convert(Integer(- v)) / m_nt_traits->convert(Integer(two * s));
      n_ys = 1;
    }
    else {
      ys_end = m_nt_traits->solve_quadratic_equation(Integer(four*r*s*s - s*t*t),
                                                     Integer(four*r*s*v - two*s*t*u),
                                                     Integer((r*v*v - t*u*v) + (t*t*w)),
                                                     ys);

      n_ys = static_cast<int>(ys_end - ys);
    }

    // Pair the x and y coordinates and obtain the vertical tangency points.
    int n(0);

    for (int i = 0; i < n_xs; ++i) {
      if (n_ys == 1) {
        ps[n++] = Point_2(xs[i], ys[0]);
      }
      else {
        for (int j = 0; j < n_ys; ++j) {
          if (CGAL::compare(m_nt_traits->convert(Integer(two*s)) * ys[j],
                            -(m_nt_traits->convert(t) * xs[i] +
                              m_nt_traits->convert(v))) == EQUAL)
          {
            ps[n++] = Point_2(xs[i], ys[j]);
            break;
          }
        }
      }
    }

    CGAL_assertion(n <= 2);
    return n;
  }

  /*! Calculate the vertical tangency points of the arc.
   * \param vpts The vertical tangency points.
   * \pre The vpts vector should be allocated at the size of 2.
   * \return The number of vertical tangency points.
   */
  size_t vertical_tangency_points(const Curve_2& cv, Alg_point_2* vpts) const {
    // No vertical tangency points for line segments:
    if (cv.orientation() == COLLINEAR) return 0;

    // Calculate the vertical tangency points of the supporting conic.
    Alg_point_2 ps[2];
    auto n = conic_vertical_tangency_points(cv, ps);
    // Return only the points that are contained in the arc interior.
    size_t m(0);
    for (int i = 0; i < n; ++i) {
      if (cv.is_full_conic() || is_strictly_between_endpoints(cv, ps[i]))
        vpts[m++] = ps[i];
    }

    // Return the number of vertical tangency points found.
    CGAL_assertion(m <= 2);
    return m;
  }

  /*! Find the horizontal tangency points of the undelying conic.
   * \param ps The output points of horizontal tangency.
   *           This area must be allocated at the size of 2.
   * \return The number of horizontal tangency points.
   */
  size_t conic_horizontal_tangency_points(const Curve_2& cv, Alg_point_2* ps)
    const {
    const Integer zero(0);

    // In case the base conic is of degree 1 (and not 2), the arc has no
    // vertical tangency points.
    if (CGAL::sign(cv.r()) == ZERO) return 0;

    // We are interested in the y coordinates were the quadratic equation:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    // has a single solution (obviously if r = 0, there are no such points).
    // We therefore demand that the discriminant of this equation is zero:
    //  (t*y + u)^2 - 4*r*(s*y^2 + v*y + w) = 0
    const Integer two(2);
    const Integer four(4);
    Algebraic ys[2];

    auto r = cv.r();
    auto s = cv.s();
    auto t = cv.t();
    auto u = cv.u();
    auto v = cv.v();
    auto w = cv.w();
    Algebraic* ys_end = m_nt_traits->template
      solve_quadratic_equation<Integer>(t*t - four*r*s,
                                        two*t*u - four*r*v,
                                        u*u - four*r*w,
                                        ys);
    auto n = static_cast<int>(ys_end - ys);

    // Compute the x coordinates and construct the horizontal tangency points.
    for (int i = 0; i < n; ++i) {
      // Having computed y, x is the single solution to the quadratic equation
      // above, and since its discriminant is 0, x is simply given by:
      Algebraic x = -(m_nt_traits->convert(t)*ys[i] + m_nt_traits->convert(u)) /
        m_nt_traits->convert(Integer(two*r));
      ps[i] = Point_2(x, ys[i]);
    }

    CGAL_assertion(n <= 2);
    return n;
  }

  /*! Calculate the horizontal tangency points of the arc.
   * \param hpts The horizontal tangency points.
   * \pre The hpts vector should be allocated at the size of 2.
   * \return The number of horizontal tangency points.
   */
  int horizontal_tangency_points(const Curve_2& cv, Alg_point_2* hpts) const {
    // No horizontal tangency points for line segments:
    if (cv.orientation() == COLLINEAR) return 0;

    // Calculate the horizontal tangency points of the conic.
    Alg_point_2 ps[2];
    std::size_t n = conic_horizontal_tangency_points(cv, ps);

    // Return only the points that are contained in the arc interior.
    int m = 0;

    for (std::size_t i = 0; i < n; ++i) {
      if (cv.is_full_conic() || is_strictly_between_endpoints(cv, ps[i]))
        hpts[m++] = ps[i];
    }

    // Return the number of horizontal tangency points found.
    CGAL_assertion(m <= 2);
    return m;
  }

  /*! Apply the inverse of the rotation given by the sin and cosine of the
   * rotation angle to the given conic arc.
   */
  void inverse_conic(const X_monotone_curve_2& xcv,
                     double cost, double sint,
                     double& r_m, double& s_m, double& t_m,
                     double& u_m, double& v_m, double& w_m) const {
    auto r = CGAL::to_double(xcv.r());
    auto s = CGAL::to_double(xcv.s());
    auto t = CGAL::to_double(xcv.t());
    auto u = CGAL::to_double(xcv.u());
    auto v = CGAL::to_double(xcv.v());
    auto w = CGAL::to_double(xcv.w());

    r_m = r * cost*cost + t*cost*sint + s*sint*sint;
    t_m = 0;
    s_m = r * sint*sint - t*cost*sint + s*cost*cost;
    u_m = u*cost + v*sint;
    v_m = - u*sint + v*cost;
    w_m = w;
  }

  /*! Obtain (i) the rotation that yields the given conic arc when applied
   * to the canonical arc, and (ii) the canonical arc.
   * \param[in] arc the given arc
   * \param[out] r_m the coefficients of the canonical conic.
   * \param[out] s_m
   * \param[out] t_m
   * \param[out] u_m
   * \param[out] v_m
   * \param[out] w_m
   * \param[out] cost the cosine of the rotation angle.
   * \param[out] sint the sine of the rotation angle.
   */
  void canonical_conic(const X_monotone_curve_2& xcv,
                       double& r_m, double& s_m, double& t_m,
                       double& u_m, double& v_m, double& w_m,
                       double& cost, double& sint) const {
    auto r = CGAL::to_double(xcv.r());
    auto s = CGAL::to_double(xcv.s());
    auto t = CGAL::to_double(xcv.t());
    // auto u = CGAL::to_double(xcv.u());
    // auto v = CGAL::to_double(xcv.v());
    // auto w = CGAL::to_double(xcv.w());
    // std::cout << r << "," << s << "," << t << ","
    //           << u << "," << v << "," << w << std::endl;

    // Compute the cos and sin of the rotation angle
    // This eliminates the t coefficinet (which multiplies x·y).

    if (r != s) {
      auto theta = atan2(t, r-s) * 0.5;
      cost = std::cos(theta);
      sint = std::sin(theta);
    }
    else if (r != 0) {
      cost = 1.0;
      sint = 0.0;
    }
    else {
      cost = std::sqrt(0.5);
      sint = cost;
    }

    inverse_conic(xcv, cost, sint, r_m, s_m, t_m, u_m, v_m, w_m);
  }

  /*! Inverse transform a point. In particular, inversly rotate the point
   * (`x`,`y`) by an angle, the sine and cosine of which are `sint` and
   * `cost`, respectively, and translate by (`-cx`,`-cy`).
   */
  void inverse_transform_point(double x, double y,
                               double cost, double sint,
                               double cx, double cy,
                               double& xc, double& yc) const {
    xc = x*cost + y*sint - cx;
    yc = -x*sint + y*cost - cy;
  }

  /*! Handle parabolas.
   * The arc-length closed form can be found here:
   * https://www.vcalc.com/wiki/vCalc/Parabola+-+arc+length
   */
  void approximate_parabola(const X_monotone_curve_2& xcv,
                            double& r_m, double& t_m, double& s_m,
                            double& u_m, double& v_m, double& w_m,
                            double& cost, double& sint,
                            double& xs_t, double& ys_t,
                            double& xt_t, double& yt_t,
                            double& a, double& ts, double& tt,
                            double& cx, double& cy,
                            bool l2r = true)
    const {
    auto min_vertex = construct_min_vertex_2_object();
    auto max_vertex = construct_max_vertex_2_object();
    const auto& src = (l2r) ? min_vertex(xcv) : max_vertex(xcv);
    const auto& trg = (l2r) ? max_vertex(xcv) : min_vertex(xcv);
    auto xs = CGAL::to_double(src.x());
    auto ys = CGAL::to_double(src.y());
    auto xt = CGAL::to_double(trg.x());
    auto yt = CGAL::to_double(trg.y());

    canonical_conic(xcv, r_m, s_m, t_m, u_m, v_m, w_m, cost, sint);
    // std::cout << r_m << "," << s_m << "," << t_m << ","
    //           << u_m << "," << v_m << "," << w_m << std::endl;
    // std::cout << "sint, cost: " << sint << "," << cost << std::endl;

    /* If the axis of the parabola is the 𝑌-axis, shift
     * the parabola by 90, essentially, converting the
     * parabola to one the axis of which is the 𝑋-axis, as the
     * remaining code assume that the 𝑋-axis is the parabola axis.
     *
     * We need to test whether s_m vanished; however, because of limited
     * precision, s_m can become very small. Therefore, instead for comparing
     * s_m with zero we compare its absolute value with the absolute value of
     * r_m, which is expected to be larger.
     */
    if (std::abs(s_m) < std::abs(r_m)) {
      // Shift phase:
      auto tmp(cost);
      cost = sint;
      sint = -tmp;

      // Recompute:
      inverse_conic(xcv, cost, sint, r_m, s_m, t_m, u_m, v_m, w_m);
    }
    // std::cout << r_m << "," << s_m << "," << t_m << ","
    //           << u_m << "," << v_m << "," << w_m << std::endl;

    // Compute the center of the inversly rotated parabola:
    // double cx_m = -u_m / (2*r_m);
    // double cy_m = (u_m*u_m - 4*r_m*w_m) / (4*r_m*v_m);
    auto cx_m = (v_m*v_m - 4*s_m*w_m) / (4*s_m*u_m);
    auto cy_m = -v_m / (2*s_m);
    // std::cout << "cx_m, cy_m: " << cx_m << "," << cy_m <<  std::endl;

    // Inverse transform the source and target
    inverse_transform_point(xs, ys, cost, sint, cx_m, cy_m, xs_t, ys_t);
    inverse_transform_point(xt, yt, cost, sint, cx_m, cy_m, xt_t, yt_t);

    a = -u_m / (4.0*s_m);
    ts = ys_t / (2.0*a);
    tt = yt_t / (2.0*a);

    // Compute the center (cx,cy) of the ellipse, rotating back:
    cx = cx_m*cost - cy_m*sint;
    cy = cx_m*sint + cy_m*cost;
    // std::cout << "center: " << cx << "," << cy << std::endl;
  }

  /*! Handle ellipses.
   */
  void approximate_ellipse(const X_monotone_curve_2& xcv,
                           double& r_m, double& t_m, double& s_m,
                           double& u_m, double& v_m, double& w_m,
                           double& cost, double& sint,
                           double& xs_t, double& ys_t, double& ts,
                           double& xt_t, double& yt_t, double& tt,
                           double& a, double& b, double& cx, double& cy,
                           bool l2r = true)
    const {
    auto min_vertex = construct_min_vertex_2_object();
    auto max_vertex = construct_max_vertex_2_object();
    const auto& src = (l2r) ? min_vertex(xcv) : max_vertex(xcv);
    const auto& trg = (l2r) ? max_vertex(xcv) : min_vertex(xcv);
    auto xs = CGAL::to_double(src.x());
    auto ys = CGAL::to_double(src.y());
    auto xt = CGAL::to_double(trg.x());
    auto yt = CGAL::to_double(trg.y());

    canonical_conic(xcv, r_m, s_m, t_m, u_m, v_m, w_m, cost, sint);
    // std::cout << r_m << "," << s_m << "," << t_m << ","
    //           << u_m << "," << v_m << "," << w_m << std::endl;
    // std::cout << "sint, cost: " << sint << "," << cost << std::endl;

    // Compute the radi of the ellipse:
    auto numerator = -4*w_m*r_m*s_m + s_m*u_m*u_m + r_m*v_m*v_m;
    auto a_sqr = numerator / (4*r_m*r_m*s_m);
    auto b_sqr = numerator / (4*r_m*s_m*s_m);
    if (a_sqr < b_sqr) {
      // Shift phase:
      auto tmp(cost);
      cost = sint;
      sint = -tmp;

      // Recompute:
      inverse_conic(xcv, cost, sint, r_m, s_m, t_m, u_m, v_m, w_m);
      numerator = -4*w_m*r_m*s_m + s_m*u_m*u_m + r_m*v_m*v_m;
      a_sqr = numerator / (4*r_m*r_m*s_m);
      b_sqr = numerator / (4*r_m*s_m*s_m);
    }

    a = std::sqrt(a_sqr);
    b = std::sqrt(b_sqr);

    // Compute the center of the inversly rotated ellipse:
    auto cx_m = -u_m / (2*r_m);
    auto cy_m = -v_m / (2*s_m);

    // Compute the center (cx,cy) of the ellipse, rotating back:
    cx = cx_m*cost - cy_m*sint;
    cy = cx_m*sint + cy_m*cost;
    // std::cout << "center: " << cx << "," << cy << std::endl;

    // Inverse transform the source and target
    inverse_transform_point(xs, ys, cost, sint, cx_m, cy_m, xs_t, ys_t);
    inverse_transform_point(xt, yt, cost, sint, cx_m, cy_m, xt_t, yt_t);
    // std::cout << "xs_t,ys_t: " << xs_t << "," << ys_t << std::endl;
    // std::cout << "xt_t,yt_t: " << xt_t << "," << yt_t << std::endl;

    // Compute the parameters ts and tt such that
    // source == (x(ts),y(ts)), and
    // target == (x(tt),y(tt))
    ts = std::atan2(a*ys_t, b*xs_t);
    if (ts < 0) ts += 2*CGAL_PI;
    tt = std::atan2(a*yt_t, b*xt_t);
    if (tt < 0) tt += 2*CGAL_PI;
    auto orient(xcv.orientation());
    if (xcv.source() != src) orient = CGAL::opposite(orient);
    if (orient == COUNTERCLOCKWISE) {
      if (tt < ts) tt += 2*CGAL_PI;
    }
    else {
      if (ts < tt) ts += 2*CGAL_PI;
    }
    // std::cout << "ts,tt: " << ts << "," << tt << std::endl;
  }

  /*! Handle hyperbolas.
   */
  void approximate_hyperbola(const X_monotone_curve_2& xcv,
                             double& r_m, double& t_m, double& s_m,
                             double& u_m, double& v_m, double& w_m,
                             double& cost, double& sint,
                             double& xs_t, double& ys_t, double& ts,
                             double& xt_t, double& yt_t, double& tt,
                             double& a, double& b, double& cx, double& cy,
                             bool l2r = true)
    const {
    auto min_vertex = construct_min_vertex_2_object();
    auto max_vertex = construct_max_vertex_2_object();
    const auto& src = (l2r) ? min_vertex(xcv) : max_vertex(xcv);
    const auto& trg = (l2r) ? max_vertex(xcv) : min_vertex(xcv);
    auto xs = CGAL::to_double(src.x());
    auto ys = CGAL::to_double(src.y());
    auto xt = CGAL::to_double(trg.x());
    auto yt = CGAL::to_double(trg.y());
    // std::cout << "curve: (" << xs << "," << ys
    //           << ") => (" << xt << "," << yt << ")"
    //           << std::endl;

    // If the hyperbola conjugate axis is the Y-axis, add
    canonical_conic(xcv, r_m, s_m, t_m, u_m, v_m, w_m, cost, sint);
    // std::cout << r_m << "," << s_m << "," << t_m << ","
    //           << u_m << "," << v_m << "," << w_m << std::endl;
    // std::cout << "sint, cost: " << sint << "," << cost << std::endl;

    auto numerator = -4*w_m*r_m*s_m + s_m*u_m*u_m + r_m*v_m*v_m;
    auto a_sqr = numerator / (4*r_m*r_m*s_m);
    auto b_sqr = -numerator / (4*r_m*s_m*s_m);

    /* If the conjugate axis of the canonical hyperbula is the 𝑌-axis, shift
     * the canonical hyperbula by 90, essentially, converting the canonical
     * hyperbula to one the conjugate axis of which is the 𝑋-axis, as the
     * remaining code assume that the conjugate axis is the 𝑋-axis.
     * Here,
     * 1. a_sqr = -a_sqr, b_sqr = -b_sqr, and
     * 2. x(t),y(t) = a*sinh(t), b*cosh(t)
     */
    if (a_sqr < 0) {
      // Shift phase:
      auto tmp(cost);
      cost = sint;
      sint = -tmp;

      // Recompute:
      inverse_conic(xcv, cost, sint, r_m, s_m, t_m, u_m, v_m, w_m);
      numerator = -4*w_m*r_m*s_m + s_m*u_m*u_m + r_m*v_m*v_m;
      a_sqr = numerator / (4*r_m*r_m*s_m);
      b_sqr = -numerator / (4*r_m*s_m*s_m);
    }
    // std::cout << "sint, cost: " << sint << "," << cost << std::endl;

    // Compute the center of the inversly rotated ellipse:
    auto cx_m = -u_m / (2*r_m);
    auto cy_m = -v_m / (2*s_m);

    // Inverse transform the source and target
    inverse_transform_point(xs, ys, cost, sint, cx_m, cy_m, xs_t, ys_t);
    inverse_transform_point(xt, yt, cost, sint, cx_m, cy_m, xt_t, yt_t);
    // std::cout << "xs_t,ys_t: " << xs_t << "," << ys_t << std::endl;
    // std::cout << "xt_t,yt_t: " << xt_t << "," << yt_t << std::endl;

    // Compute the center (cx,cy) of the hyperbola, rotating back:
    cx = cx_m*cost - cy_m*sint;
    cy = cx_m*sint + cy_m*cost;
    // std::cout << "center: " << cx << "," << cy << std::endl;

    a = std::sqrt(a_sqr);
    b = std::sqrt(b_sqr);

    // We use the parametric representation x(t),y(t) = a*cosh(t), b*sinh(t)
    // cosh(t) = (e^t + e^(-t))/2
    // sinh(t) = (e^t - e^(-t))/2
    // Compute the parameters ts and tt such that
    // source == (x(ts),y(ts)), and
    // target == (x(tt),y(tt))
    // Compute the radi of the hyperbola:
    ts = std::asinh(ys_t/b);
    tt = std::asinh(yt_t/b);
    assert(std::signbit(xs_t) == std::signbit(xt_t));

    if (std::signbit(xs_t)) a = -a;
  }
};

#include <CGAL/enable_warnings.h>

} //namespace CGAL

#endif
