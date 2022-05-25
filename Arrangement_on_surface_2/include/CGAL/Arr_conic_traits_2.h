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

#ifndef CGAL_ARR_CONIC_TRAITS_2_H
#define CGAL_ARR_CONIC_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The conic traits-class for the arrangement package.
 */

#include <fstream>
#include <atomic>
#include <memory>
#include <cmath>

#include <CGAL/Cartesian.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_geometry_traits/Conic_arc_2.h>
#include <CGAL/Arr_geometry_traits/Conic_x_monotone_arc_2.h>
#include <CGAL/Arr_geometry_traits/Conic_point_2.h>

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

  typedef typename Nt_traits::Integer     Integer;

  typedef Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>  Self;

  // Category tags:
  typedef Tag_true                        Has_left_category;
  typedef Tag_true                        Has_merge_category;
  typedef Tag_false                       Has_do_intersect_category;
  //typedef boost::true_type                Has_line_segment_constructor;

  typedef Arr_oblivious_side_tag          Left_side_category;
  typedef Arr_oblivious_side_tag          Bottom_side_category;
  typedef Arr_oblivious_side_tag          Top_side_category;
  typedef Arr_oblivious_side_tag          Right_side_category;

  // Traits objects:
  typedef _Conic_arc_2<Rat_kernel, Alg_kernel, Nt_traits> Curve_2;
  typedef _Conic_x_monotone_arc_2<Curve_2>                X_monotone_curve_2;
  typedef _Conic_point_2<Alg_kernel>                      Point_2;
  typedef unsigned int                                    Multiplicity;

private:
  // Type definition for the intersection points mapping.
  typedef typename X_monotone_curve_2::Conic_id           Conic_id;
  typedef typename X_monotone_curve_2::Intersection_point Intersection_point;
  typedef typename X_monotone_curve_2::Intersection_map   Intersection_map;

  typedef std::shared_ptr<Rat_kernel>           Shared_rat_kernel;
  typedef std::shared_ptr<Alg_kernel>           Shared_alg_kernel;
  typedef std::shared_ptr<Nt_traits>            Shared_nt_Traits;

  const Shared_rat_kernel m_rat_kernel;
  const Shared_alg_kernel m_alg_kernel;
  const Shared_nt_Traits m_nt_traits;

  mutable Intersection_map inter_map;   // Mapping conic pairs to their
                                        // intersection points.

public:
  /*! Default constructor.
   */
  Arr_conic_traits_2() {}

  /*! Construct from resources.
   */
  Arr_conic_traits_2(Shared_rat_kernel rat_kernel,
                     Shared_alg_kernel alg_kernel,
                     Shared_nt_Traits nt_traits) :
    m_rat_kernel(rat_kernel),
    m_alg_kernel(alg_kernel),
    m_nt_traits(nt_traits)
  {}

  /*! Obtain the next conic index. */
  static unsigned int get_index() {
#ifdef CGAL_NO_ATOMIC
    static unsigned int index;
#else
    static std::atomic<unsigned int> index;
#endif
    return (++index);
  }

  /// \name Basic functor definitions.
  //@{

  class Compare_x_2 {
  public:
    /*! Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) < x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      Alg_kernel   ker;
      return (ker.compare_x_2_object()(p1, p2));
    }
  };

  /*! Obtain a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object () const { return Compare_x_2(); }

  class Compare_xy_2 {
  public:
    /*! Compares two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2);
     *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL if the two points are equal.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      Alg_kernel ker;
      return ker.compare_xy_2_object()(p1, p2);
    }
  };

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(); }

  class Construct_min_vertex_2 {
  public:
    /*! Obtain the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2 & cv) const
    { return cv.left(); }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  class Construct_max_vertex_2 {
  public:
    /*! Obtain the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The right endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2 & cv) const
    { return cv.right(); }
  };

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  class Is_vertical_2 {
  public:
    /*! Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv) const
    { return cv.is_vertical(); }
  };

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(); }

  class Compare_y_at_x_2 {
  public:
    /*! Return the location of the given point with respect to the input curve.
     * \param cv The curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & cv) const
    {
      Alg_kernel ker;

      if (cv.is_vertical()) {
        // A special treatment for vertical segments:
        // In case p has the same x c-ordinate of the vertical segment, compare
        // it to the segment endpoints to determine its position.
        Comparison_result res1 = ker.compare_y_2_object()(p, cv.left());
        Comparison_result res2 = ker.compare_y_2_object()(p, cv.right());

        if (res1 == res2) return res1;
        else return EQUAL;
      }

      // Check whether the point is exactly on the curve.
      if (cv.contains_point(p)) return EQUAL;

      // Obtain a point q on the x-monotone arc with the same x coordinate as p.
      Comparison_result x_res;
      Point_2 q;

      if ((x_res = ker.compare_x_2_object()(p, cv.left())) == EQUAL) {
        q = cv.left();
      }
      else {
        CGAL_precondition (x_res != SMALLER);

        if ((x_res = ker.compare_x_2_object()(p, cv.right())) == EQUAL) {
          q = cv.right();
        }
        else {
          CGAL_precondition(x_res != LARGER);

          q = cv.point_at_x (p);
        }
      }

      // Compare p with the a point of the curve with the same x coordinate.
      return (ker.compare_y_2_object()(p, q));
    }
  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(); }

  class Compare_y_at_x_left_2 {
  public:
    /*! Compares the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& p) const
    {
      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition(cv1.contains_point(p) &&
                        cv2.contains_point(p));

      CGAL_precondition_code(Alg_kernel ker;);
      CGAL_precondition(ker.compare_xy_2_object()(p, cv1.left()) == LARGER &&
                        ker.compare_xy_2_object()(p, cv2.left()) == LARGER);

      // If one of the curves is vertical, it is below the other one.
      if (cv1.is_vertical()) {
        // Check whether both are vertical:
        if (cv2.is_vertical()) return EQUAL;
        else return SMALLER;
      }
      else if (cv2.is_vertical()) return LARGER;

      // Compare the two curves immediately to the left of p:
      return cv1.compare_to_left(cv2, p);
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(); }

  class Compare_y_at_x_right_2 {
  public:
    /*! Compares the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& p) const
    {
      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition(cv1.contains_point(p) && cv2.contains_point(p));
      CGAL_precondition_code(Alg_kernel ker;);
      CGAL_precondition(ker.compare_xy_2_object()(p, cv1.right()) == SMALLER &&
                        ker.compare_xy_2_object()(p, cv2.right()) == SMALLER);

      // If one of the curves is vertical, it is above the other one.
      if (cv1.is_vertical()) {
        // Check whether both are vertical:
        if (cv2.is_vertical()) return EQUAL;
        else return LARGER;
      }
      else if (cv2.is_vertical()) return SMALLER;

      // Compare the two curves immediately to the right of p:
      return cv1.compare_to_right(cv2, p);
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(); }

  class Equal_2 {
  public:
    /*! Check whether the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      if (&cv1 == &cv2) return true;
      return cv1.equals(cv2);
    }

    /*! Check whether the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const {
      if (&p1 == &p2) return (true);

      Alg_kernel ker;
      return(ker.compare_xy_2_object()(p1, p2) == EQUAL);
    }
  };

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object() const { return Equal_2(); }
  //@}

  /// \name Intersections, subdivisions, and mergings
  //@{

  /*! \class Make_x_monotone_2
   * A functor for subdividing curves into x-monotone curves.
   */
  class Make_x_monotone_2 {
    typedef Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>       Self;

  public:
    /*! Subdivide a given conic curve (or conic arc) into x-monotone subcurves
     * and insert them to a given output iterator.
     * \param cv the curve.
     * \param oi the output iterator for the result. Its dereference type is a
     *           variant that wraps a \c Point_2 or an \c X_monotone_curve_2
     *           objects.
     * \return the past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {
      typedef boost::variant<Point_2, X_monotone_curve_2>
        Make_x_monotone_result;

      // Increment the serial number of the curve cv, which will serve as its
      // unique identifier.
      auto index = Self::get_index();
      Conic_id conic_id(index);

      // Find the points of vertical tangency to cv and act accordingly.
      typename Curve_2::Point_2 vtan_ps[2];
      int n_vtan_ps;

      n_vtan_ps = cv.vertical_tangency_points(vtan_ps);

      if (n_vtan_ps == 0) {
        // In case the given curve is already x-monotone:
        *oi++ = Make_x_monotone_result(X_monotone_curve_2(cv, conic_id));
        return oi;
      }

      // Split the conic arc into x-monotone sub-curves.
      if (cv.is_full_conic()) {
        // Make sure we have two vertical tangency points.
        CGAL_assertion(n_vtan_ps == 2);

        // In case the curve is a full conic, split it into two x-monotone
        // arcs, one going from ps[0] to ps[1], and the other from ps[1] to
        // ps[0].
        *oi++ = Make_x_monotone_result(X_monotone_curve_2(cv, vtan_ps[0],
                                                          vtan_ps[1],
                                                          conic_id));
        *oi++ = Make_x_monotone_result(X_monotone_curve_2(cv, vtan_ps[1],
                                                          vtan_ps[0],
                                                          conic_id));
      }
      else {
        if (n_vtan_ps == 1) {
          // Split the arc into two x-monotone sub-curves: one going from the
          // arc source to ps[0], and the other from ps[0] to the target.
          *oi++ = Make_x_monotone_result(X_monotone_curve_2(cv, cv.source(),
                                                            vtan_ps[0],
                                                            conic_id));
          *oi++ = Make_x_monotone_result(X_monotone_curve_2(cv, vtan_ps[0],
                                                            cv.target(),
                                                            conic_id));
        }
        else {
          CGAL_assertion(n_vtan_ps == 2);

          // Identify the first point we encounter when going from cv's source
          // to its target, and the second point we encounter. Note that the
          // two endpoints must both be below the line connecting the two
          // tangnecy points (or both lies above it).
          int ind_first = 0;
          int ind_second = 1;
          Alg_kernel ker;
          typename Alg_kernel::Line_2 line =
            ker.construct_line_2_object()(vtan_ps[0], vtan_ps[1]);
          const Comparison_result start_pos =
            ker.compare_y_at_x_2_object()(cv.source(), line);
          const Comparison_result order_vpts =
            ker.compare_x_2_object()(vtan_ps[0], vtan_ps[1]);

          CGAL_assertion(start_pos != EQUAL &&
                         ker.compare_y_at_x_2_object()(cv.target(),
                                                       line) == start_pos);
          CGAL_assertion(order_vpts != EQUAL);

          if (((cv.orientation() == COUNTERCLOCKWISE) &&
               (start_pos == order_vpts)) ||
              ((cv.orientation() == CLOCKWISE) && (start_pos != order_vpts)))
          {
            ind_first = 1;
            ind_second = 0;
          }

          // Split the arc into three x-monotone sub-curves.
          *oi++ = Make_x_monotone_result(X_monotone_curve_2(cv, cv.source(),
                                                            vtan_ps[ind_first],
                                                            conic_id));

          *oi++ = Make_x_monotone_result(X_monotone_curve_2(cv,
                                                            vtan_ps[ind_first],
                                                            vtan_ps[ind_second],
                                                            conic_id));

          *oi++ = Make_x_monotone_result(X_monotone_curve_2(cv,
                                                            vtan_ps[ind_second],
                                                            cv.target(),
                                                            conic_id));
        }
      }

      return oi;
    }
  };

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(); }

  class Split_2 {
  public:
    /*! Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint).
     * \param c2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& cv, const Point_2 & p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    { cv.split(p, c1, c2); }
  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object() const { return Split_2(); }

  class Intersect_2 {
  private:
    Intersection_map& _inter_map;       // The map of intersection points.

  public:
    /*! Constructor. */
    Intersect_2(Intersection_map& map) : _inter_map(map) {}

    /*! Find the intersections of the two given curves and insert them to the
     * given output iterator. As two segments may itersect only once, only a
     * single will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    { return cv1.intersect(cv2, _inter_map, oi); }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const { return (Intersect_2(inter_map)); }

  class Are_mergeable_2 {
  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    { return cv1.can_merge_with(cv2); }
  };

  /*! Obtain an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(); }

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2 {
  protected:
    typedef Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>       Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

  public:
    /*! Merge two given x-monotone curves into a single curve (segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable.
     */
    void operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const
    {
      CGAL_precondition(m_traits->are_mergeable_2_object()(cv2, cv1));

      c = cv1;
      c.merge (cv2);
    }
  };

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object() const { return Merge_2(this); }

  //@}

  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  using Approximate_number_type = double;
  using Approximate_kernel = CGAL::Cartesian<Approximate_number_type>;
  using Approximate_point_2 = Approximate_kernel::Point_2;

  class Approximate_2 {
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
    { return std::make_pair(operator()(p, 0), operator()(p, 1)); }

    /*! Obtain an approximation of an x-monotone curve.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& arc, size_t size,
                              OutputIterator oi) const {
      auto xs = CGAL::to_double(arc.source().x());
      auto ys = CGAL::to_double(arc.source().y());
      auto xt = CGAL::to_double(arc.target().x());
      auto yt = CGAL::to_double(arc.target().y());
      if (arc.orientation() == COLLINEAR) {
        *oi++ = Approximate_point_2(xs, ys);
        *oi++ = Approximate_point_2(xt, yt);
        return oi;
      }
      auto r = CGAL::to_double(arc.r());
      auto s = CGAL::to_double(arc.s());
      auto t = CGAL::to_double(arc.t());
      auto u = CGAL::to_double(arc.u());
      auto v = CGAL::to_double(arc.v());
      auto w = CGAL::to_double(arc.w());
      std::cout << r << "," << s << "," << t << "," << u << "," << v << "," << w
                << std::endl;
      {
        // Compute the cos and sin of the rotation angle
        // In case of a circle, cost == 1 and sint = 0
        double cost(1), sint(0);

        if (r != s) {
          auto tan_2t = t / (r - s);
          auto cos_2t = std::sqrt(1 / (tan_2t*tan_2t + 1));
          cost = std::sqrt((1 + cos_2t) / 2);
          sint = std::sqrt((1 - cos_2t) / 2);
        }
        std::cout << "sint, cost: " << sint << "," << cost << std::endl;

        // Compute the coefficients of the unrotated ellipse
        auto r_m = r * cost*cost + t*cost*sint + s*sint*sint;
        auto t_m = 0;
        auto s_m = r * sint*sint - t*cost*sint + s*cost*cost;
        auto u_m = u*cost + v*sint;
        auto v_m = - u*sint + v*cost;
        auto w_m = w;

        std::cout << r_m << "," << s_m << "," << t_m << ","
                  << u_m << "," << v_m << "," << w_m << std::endl;

        // Compute the center of the inversly rotated ellipse:
        auto cx_m = -u_m / (2*r_m);
        auto cy_m = -v_m / (2*s_m);

        // Compute the radi of the ellipse:
        auto numerator = -4*w_m*r_m*s_m + s_m*u_m*u_m + r_m*v_m*v_m;
        auto a = std::sqrt(numerator / (4*r_m*r_m*s_m));
        auto b = std::sqrt(numerator / (4*r_m*s_m*s_m));
        std::cout << "a, b: " << a << "," << b << std::endl;

        // Compute the center (cx,cy) of the ellipse, rotating back:
        auto cx = cx_m*cost - cy_m*sint;
        auto cy = cx_m*sint + cy_m*cost;
        std::cout << "center: " << cx << "," << cy << std::endl;

        // Compute the parameters ts and tt such that
        // source == (x(ts),y(ts)), and
        // target == (x(tt),y(tt))
        auto xds = xs - cx;
        auto yds = ys - cy;
        auto ts = std::atan2(a*(cost*yds - sint*xds),b*(sint*yds + cost*xds));
        auto xdt = xt - cx;
        auto ydt = yt - cy;
        auto tt = std::atan2(a*(cost*ydt - sint*xdt),b*(sint*ydt + cost*xdt));

        auto delta = std::abs(tt - ts) / (size-1);
        double t((arc.orientation() == COUNTERCLOCKWISE) ? ts : tt);

        *oi++ = (arc.orientation() == COUNTERCLOCKWISE) ?
          Approximate_point_2(xs, ys) : Approximate_point_2(xt, yt);
        t += delta;
        for (size_t i = 1; i < size-1; ++i) {
          auto x = a*std::cos(t)*cost - b*std::sin(t)*sint + cx;
          auto y = a*std::cos(t)*sint + b*std::sin(t)*cost + cy;
          std::cout << "t, (x, y): " << t << ", (" << x << "," << y << ")"
                    << std::endl;
          *oi++ = Approximate_point_2(x, y);
          t += delta;
        }
        *oi++ = (arc.orientation() == COUNTERCLOCKWISE) ?
          Approximate_point_2(xt, yt) : Approximate_point_2(xs, ys);
      }
      return oi;
    }
  };

  /*! Obtain an Approximate_2 functor object. */
  Approximate_2 approximate_2_object() const { return Approximate_2(); }

  //! Functor
  class Construct_x_monotone_curve_2 {
  public:
    /*! Return an x-monotone curve connecting the two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q must not be the same.
     * \return A segment connecting p and q.
     */
    X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q) const
    { return (X_monotone_curve_2(p, q)); }
  };

  /*! Obtain a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object () const
  { return Construct_x_monotone_curve_2(); }

  //! Constructor of conic arcs
  class Construct_curve_2 {
  protected:
    typedef Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>  Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits.
     */
    Construct_curve_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>;

  public:
    /*! Construct an empty curve.
     */
    Curve_2 operator()() const { return Curve_2(); }

    /*! Construct a conic arc which is the full conic:
     *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
     * \pre The conic C must be an ellipse (so 4rs - t^2 > 0).
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
      set_full(rat_coeffs, true);
    }

    /*! Construct a conic arc that lies on the conic:
     *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
     * \param orient The orientation of the arc (clockwise or counterclockwise).
     * \param source The source point.
     * \param target The target point.
     * \pre The source and the target must be on the conic boundary and must
     * not be the same.
     */
    Curve_2 operator()(const Rational& r, const Rational& s, const Rational& t,
                       const Rational& u, const Rational& v, const Rational& w,
                       const Orientation& orient,
                       const Point_2& source, const Point_2& target) const {
      // Make sure that the source and the taget are not the same.
      const auto alg_kernel = m_traits.m_alg_kernel;
      CGAL_precondition(alg_kernel->compare_xy_2_object()(source, target) !=
                        EQUAL);
      // Set the arc properties (no need to compute the orientation).
      Rational rat_coeffs[6] = {r, s, t, u, v, w};
      Curve_2 arc;
      arc.set_orientation(orient);
      arc.set_endpoints(source, target);
      set(arc, rat_coeffs);
      return arc;
    }

    /*! Construct a conic arc that is a circular arc from given three points.
     * \param p1 The arc source.
     * \param p2 A point in the interior of the arc.
     * \param p3 The arc target.
     * \pre The three points must not be collinear.
     */
    Curve_2 operator()(const Rat_point_2& p1, const Rat_point_2& p2,
                       const Rat_point_2& p3) {
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
      auto source = Point_2(nt_traits->convert(x1), nt_traits->convert(y1));
      auto target = Point_2(nt_traits->convert(x3), nt_traits->convert(y3));
      arc.set_enpoints(source, target);

      // Make sure that the source and the taget are not the same.
      CGAL_precondition(alg_kernel->compare_xy_2_object()(source, target) !=
                        EQUAL);

      // Compute the lines: A1*x + B1*y + C1 = 0,
      //               and: A2*x + B2*y + C2 = 0,
      // where:
      const Rational two  = 2;

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
      auto orient_f = m_alg_kernel.orientation_2_object();
      Point_2 p_mid = Point_2(nt_traits->convert(x2), nt_traits->convert(y2));

      auto orient = (orient_f(source, p_mid, target) == LEFT_TURN) ?
        COUNTERCLOCKWISE : CLOCKWISE;
      arc.set_orientation(orient);

      // Set the arc properties (no need to compute the orientation).
      set(arc, rat_coeffs);
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
                       const Rat_point_2& p5) {
      Curve_2 arc;

      // Make sure that no three points are collinear.
      Rat_kernel ker;
      auto orient_f = m_rat_kernel.orientation_2_object();
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
      auto source = Point_2(nt_traits->convert(x1), nt_traits->convert(y1));
      auto target = Point_2(nt_traits->convert(x5), nt_traits->convert(y5));
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
      set(arc, rat_coeffs);

      // Make sure that all midpoints are strictly between the
      // source and the target.
      Point_2 mp2 =
        Point_2(nt_traits->convert(p2.x()), nt_traits->convert(p2.y()));
      Point_2 mp3 =
        Point_2(nt_traits->convert(p3.x()), m_nt_traits->convert(p3.y()));
      Point_2 mp4 =
        Point_2(nt_traits->convert(p4.x()), nt_traits->convert(p4.y()));

      if (! is_strictly_between_endpoints(arc, mp2) ||
          ! is_strictly_between_endpoints(arc, mp3) ||
          ! is_strictly_between_endpoints(arc, mp4))
      {
        arc.reset_flags();            // inavlid arc
        return arc;
      }
      return arc;
    }

    /*! Construct a conic arc that lies on a conic given by its coefficients:
     *   C: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
     * The source and the target are specified by the intersection of the
     * conic with:
     *   C_1: r_1*x^2 + s_1*y^2 + t_1*xy + u_1*x + v_1*y + w_1 = 0
     *   C_2: r_2*x^2 + s_2*y^2 + t_2*xy + u_2*x + v_2*y + w_2 = 0
     * The user must also specify the source and the target with approximated
     * coordinates. The actual intersection points that best fits the source
     * (or the target) will be selected.
     */
    Curve_2 operator()(const Rational& r, const Rational& s, const Rational& t,
                       const Rational& u, const Rational& v, const Rational& w,
                       const Orientation& orient,
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
        // Get the integer coefficients of the k'th auxiliary conic curve.
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
        // conic and the k'th auxiliary conic.
        int n_xs = compute_resultant_roots(nt_traits,
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

        int n_ys = compute_resultant_roots(nt_traits,
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
      set(arc, rat_coeffs);
      return arc;
    }

    /*! Return a curve connecting the two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q must not be the same.
     * \return A segment connecting p and q.
     */
    Curve_2 operator()(const Point_2& p, const Point_2& q) const {
      const auto alg_kernel = m_traits.m_alg_kernel;
      CGAL_precondition(alg_kernel->compare_xy_2_object()(p, q) != EQUAL);

      Curve_2 arc;
      set_orientation(Curve_2::COLLINEAR);
      arc.set_flag(Curve_2::IS_VALID);
      arc.set_source(p);
      arc.set_target(q);

      // Compose the equation of the underlying line.
      const Algebraic& x1 = p.x();
      const Algebraic& y1 = p.y();
      const Algebraic& x2 = q.x();
      const Algebraic& y2 = q.y();

      // The supporting line is A*x + B*y + C = 0, where:
      //  A = y2 - y1,    B = x1 - x2,    C = x2*y1 - x1*y2
      // We use the extra data field to store the equation of this line.
      auto* extra_data = new typename Curve_2::Extra_data;
      extra_data->a = y2 - y1;
      extra_data->b = x1 - x2;
      extra_data->c = x2*y1 - x1*y2;
      extra_data->side = ZERO;
      arc.set_extra_data(extra_data);
      return arc;
    }

    /*! Construct a conic arc from a given line segment.
     * \param seg The line segment with rational endpoints.
     */
    Curve_2 operator()(const Rat_segment_2& seg) const {
      Curve_2 arc;
      arc.set_orientation(COLLINEAR);

      // Set the source and target.
      const auto rat_kernel = m_traits.m_rat_kernel;
      Rat_point_2 source = rat_kernel->construct_vertex_2_object()(seg, 0);
      Rat_point_2 target = rat_kernel->construct_vertex_2_object()(seg, 1);
      const Rational& x1 = source.x();
      const Rational& y1 = source.y();
      const Rational& x2 = target.x();
      const Rational& y2 = target.y();

      const auto nt_traits = m_traits.m_nt_traits;
      arc.set_source(Point_2(nt_traits->convert(x1), nt_traits->convert(y1)));
      arc.set_target(Point_2(nt_traits->convert(x2), nt_traits->convert(y2)));

      // Make sure that the source and the taget are not the same.
      CGAL_precondition(rat_kernel->compare_xy_2_object()(source, target) !=
                        EQUAL);

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
      m_traits.set(arc, rat_coeffs);
      return arc;
    }

    /*! Construct a conic arc that is a full circle.
     * \param circ The circle with rational center and rational squared radius.
     */
    Curve_2 operator()(const Rat_circle_2& circ) const {
      Curve_2 arc;
      arc.set_orientation(CLOCKWISE);

      // Get the circle properties.
      const auto rat_kernel = m_traits.m_rat_kernel;
      Rat_point_2 center = rat_kernel->construct_center_2_object()(circ);
      Rational x0 = center.x();
      Rational y0 = center.y();
      Rational R_sqr = rat_kernel->compute_squared_radius_2_object()(circ);

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
      rat_coeffs[5] = x0*x0 + y0*y0 - R_sqr;

      // Set the arc to be the full conic (no need to compute the orientation).
      m_traits.set_full(arc.rat_coeffs, false);
      return arc;
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
    Comparison_result operator() (const X_monotone_curve_2& cv) const {
      if (cv.is_directed_right()) return SMALLER;
      else return LARGER;
    }
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
    typedef Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>       Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Trim_2(const Traits& traits) : m_traits(traits) {}

  public:
    friend class Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
    /*!\brief
     * Returns a trimmed version of an arc
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
                                  const Point_2& src,
                                  const Point_2& tgt)const
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
      if( (xcv.is_directed_right() && compare_x_2(src, tgt) == LARGER) ||
          (! xcv.is_directed_right() && compare_x_2(src, tgt) == SMALLER) )
        return (xcv.trim(tgt, src));
      else return (xcv.trim(src, tgt));
    }
  };

  /*! Obtain a Trim_2 functor object. */
  Trim_2 trim_2_object() const { return Trim_2(*this); }
  //@}

  /*! Set the properties of a conic arc (for the usage of the constructors).
   * \param rat_coeffs A vector of size 6, storing the rational coefficients
   *                   of x^2, y^2, xy, x, y and the free coefficient resp.
   */
  void set(Curve_2& arc, const Rational* rat_coeffs) const {
    arc.set_flag(Curve_2::IS_VALID);

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
    if (arc.orientation() == temp_conic.orientation()) {
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
    arc.set_coefficients(r, s, t, u, v, w);

    const auto& source = arc.source();
    const auto& target = arc.target();
    // Make sure both endpoint lie on the supporting conic.
    if (! is_on_supporting_conic(arc, source) ||
        ! is_on_supporting_conic(arc, target))
    {
      arc.reset_flags();            // inavlid arc
    }

    // Check whether we have a degree 2 curve.
    if ((CGAL::sign(r) != ZERO) || (CGAL::sign(s) != ZERO) ||
        (CGAL::sign(t) != ZERO))
    {
      if (arc.orientation() == COLLINEAR) {
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
          arc.reset_flags();    // inavlid arc
          return;
        }


        // We have a segment of a line pair with rational coefficients.
        // Compose the equation of the underlying line
        // (with algebraic coefficients).
        const Algebraic& x1 = source.x();
        const Algebraic& y1 = source.y();
        const Algebraic& x2 = target.x();
        const Algebraic& y2 = target.y();

        // The supporting line is A*x + B*y + C = 0, where:
        //
        //  A = y2 - y1,    B = x1 - x2,    C = x2*y1 - x1*y2
        //
        // We use the extra dat field to store the equation of this line.
        auto* extra_data = new typename Curve_2::Extra_data;
        extra_data->a = y2 - y1;
        extra_data->b = x1 - x2;
        extra_data->c = x2*y1 - x1*y2;
        extra_data->side = ZERO;
        arc.set_extra_data(extra_data);
      }
      else {
        // The sign of (4rs - t^2) detetmines the conic type:
        // - if it is possitive, the conic is an ellipse,
        // - if it is negative, the conic is a hyperbola,
        // - if it is zero, the conic is a parabola.
        CGAL::Sign sign_conic = CGAL::sign(4*r*s - t*t);

        // Build the extra hyperbolic data if necessary
        if (sign_conic == NEGATIVE) build_hyperbolic_arc_data(arc);

        if (sign_conic != POSITIVE) {
          // In case of a non-degenerate parabola or a hyperbola, make sure
          // the arc is not infinite.
          Point_2 p_mid =
            m_alg_kernel->construct_midpoint_2_object()(source, target);
          Point_2 ps[2];

          bool finite_at_x = (points_at_x(arc, p_mid, ps) > 0);
          bool finite_at_y = (points_at_y(arc, p_mid, ps) > 0);

          if (! finite_at_x && ! finite_at_y) {
            arc.reset_flags();              // inavlid arc
            return;
          }
        }
      }
    }


    arc.set_flag(Curve_2::IS_VALID);            // arc is valid
    arc.reset_flag(Curve_2::IS_FULL_CONIC);     // not a full conic
  }

  /*! Set the properties of a conic arc that is really a full curve
   * (that is, an ellipse).
   * \param rat_coeffs A vector of size 6, storing the rational coefficients
   *                   of x^2, y^2, xy, x, y and the free coefficient resp.
   * \param comp_orient Should we compute the orientation of the given curve.
   */
  void set_full(Curve_2& arc, const Rational* rat_coeffs,
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

    if (comp_orient) arc.set_orientation(temp_orient);

    Integer r, s, t, u, v, w;
    if (arc.orientation() == temp_orient) {
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

    // Make sure the conic is a non-degenerate ellipse:
    // The coefficients should satisfy (4rs - t^2) > 0.
    const bool is_ellipse = (CGAL::sign(4*r*s - t*t) == POSITIVE);
    CGAL_assertion(is_ellipse);

    // We do not have to store any extra data with the arc.

    // Mark that this arc is a full conic curve.
    if (is_ellipse) {
      arc.set_flag(Curve_2::IS_VALID);
      arc.set_flag(Curve_2::IS_FULL_CONIC);
    }
    else arc.reset_flags();            // inavlid arc
  }

  /*! Check whether the given point lies on the supporting conic of the arc.
   * \param p The query point.
   * \return true if p lies on the supporting conic; (false) otherwise.
   */
  bool is_on_supporting_conic(Curve_2& arc, const Point_2& p) const {
    // Check whether p satisfies the conic equation.
    const Algebraic r = m_nt_traits->convert(arc.r());
    const Algebraic t = m_nt_traits->convert(arc.t());
    const Algebraic u = m_nt_traits->convert(arc.u());
    const Algebraic s = m_nt_traits->convert(arc.s());
    const Algebraic v = m_nt_traits->convert(arc.v());
    const Algebraic w = m_nt_traits->convert(arc.w());

    // The point must satisfy: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0.
    const Algebraic val =
      (r*p.x() + t*p.y() + u)*p.x() + (s*p.y() + v)* p.y() + w;
    return (CGAL::sign(val) == ZERO);
  }

  /*! Check whether the given point is between the source and the target.
   * The point is assumed to be on the conic's boundary.
   * \param p The query point.
   * \return true if the point is between the two endpoints,
   *         (false) if it is not.
   */
  bool is_between_endpoints(const Curve_2& arc, const Point_2& p) const {
    CGAL_precondition(! arc.is_full_conic());

    // Check if p is one of the endpoints.
    auto eq = m_alg_kernel->equal_2_object();
    if (eq(p, arc.source()) || eq(p, arc.target())) return true;
    else return (is_strictly_between_endpoints(arc, p));
  }

  /*! Check whether the given point is strictly between the source and the
   * target (but not any of them).
   * The point is assumed to be on the conic's boundary.
   * \param p The query point.
   * \return true if the point is strictly between the two endpoints,
   *         (false) if it is not.
   */
  bool is_strictly_between_endpoints(const Curve_2& arc, const Point_2& p) const
  {
    // In case this is a full conic, any point on its boundary is between
    // its end points.
    if (arc.is_full_conic()) return true;

    // Check if we have extra data available.
    const auto* extra_data = arc.extra_data();
    if (extra_data != nullptr) {
      if (extra_data->side != ZERO) {
        // In case of a hyperbolic arc, make sure the point is located on the
        // same branch as the arc.
        if (arc.sign_of_extra_data(p.x(), p.y()) != extra_data->side)
          return false;
      }
      else {
        // In case we have a segment of a line pair, make sure that p really
        // satisfies the equation of the line.
        if (arc.sign_of_extra_data(p.x(), p.y()) != ZERO) return false;
      }
    }

    // Act according to the conic degree.
    auto orient_f = m_alg_kernel->orientation_2_object();
    const auto& source = arc.source();
    const auto& target = arc.target();
    if (arc.orientation() == COLLINEAR) {
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
    else {
      // In case of a conic of degree 2, make a decision based on the conic's
      // orientation and whether (source,p,target) is a right or a left turn.
      if (arc.orientation() == COUNTERCLOCKWISE)
        return (orient_f(source, p, target) == LEFT_TURN);
      else
        return (orient_f(source, p, target) == RIGHT_TURN);
    }
  }

  /*! Build the data for hyperbolic arc, contaning the characterization of the
   * hyperbolic branch the arc is placed on.
   */
  void build_hyperbolic_arc_data(Curve_2& arc) const {
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
    const int or_fact = (arc.orientation() == CLOCKWISE) ? -1 : 1;
    const Algebraic r = m_nt_traits->convert(or_fact * arc.r());
    const Algebraic s = m_nt_traits->convert(or_fact * arc.s());
    const Algebraic t = m_nt_traits->convert(or_fact * arc.t());
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
    const Algebraic u = m_nt_traits->convert(or_fact * arc.u());
    const Algebraic v = m_nt_traits->convert(or_fact * arc.v());
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
    auto* extra_data = new typename Curve_2::Extra_data;

    extra_data->a = cos_phi;
    extra_data->b = sin_phi;
    extra_data->c = - (cos_phi*x0 + sin_phi*y0);

    // Make sure that the two endpoints are located on the same branch
    // of the hyperbola.
    const auto& source = arc.source();
    const auto& target = arc.target();
    extra_data->side = arc.sign_of_extra_data(source.x(), source.y());

    CGAL_assertion(extra_data->side != ZERO);
    CGAL_assertion(extra_data->side ==
                   arc.sign_of_extra_data(target.x(), target.y()));
    arc.set_extra_data(extra_data);
  }

  /*! Find the x coordinates of the underlying conic at a given y coordinate.
   * \param y The y coordinate.
   * \param xs The output x coordinates.
   * \pre The vector xs must be allocated at the size of 2.
   * \return The number of x coordinates computed (either 0, 1 or 2).
   */
  int conic_get_x_coordinates(const Curve_2& arc,
                              const Algebraic& y, Algebraic* xs) const {
    // Solve the quadratic equation for a given y and find the x values:
    //  r*x^2 + (t*y + u)*x + (s*y^2 + v*y + w) = 0
    Algebraic A = m_nt_traits->convert(arc.r());
    Algebraic B =
      m_nt_traits->convert(arc.t())*y + m_nt_traits->convert(arc.u());
    Algebraic C =
      (m_nt_traits->convert(arc.s())*y + m_nt_traits->convert(arc.v()))*y +
      m_nt_traits->convert(arc.w());

    return solve_quadratic_equation(A, B, C, xs[0], xs[1]);
  }

  /*! Find the y coordinates of the underlying conic at a given x coordinate.
   * \param x The x coordinate.
   * \param ys The output y coordinates.
   * \pre The vector ys must be allocated at the size of 2.
   * \return The number of y coordinates computed (either 0, 1 or 2).
   */
  int conic_get_y_coordinates(const Curve_2& arc,
                              const Algebraic& x, Algebraic* ys) const {
    // Solve the quadratic equation for a given x and find the y values:
    //  s*y^2 + (t*x + v)*y + (r*x^2 + u*x + w) = 0
    Algebraic A = m_nt_traits->convert(arc.s());
    Algebraic B =
      m_nt_traits->convert(arc.t())*x + m_nt_traits->convert(arc.v());
    Algebraic C =
      (m_nt_traits->convert(arc.r())*x + m_nt_traits->convert(arc.u()))*x +
      m_nt_traits->convert(arc.w());

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
    Nt_traits nt_traits;
    Algebraic sqrt_disc = m_nt_traits->sqrt(disc);

    x_minus = -(B + sqrt_disc) / _2A;
    x_plus = (sqrt_disc - B) / _2A;
    return 2;
  }

  /*! Find all points on the arc with a given x-coordinate.
   * \param p A placeholder for the x-coordinate.
   * \param ps The point on the arc at x(p).
   * \pre The vector ps should be allocated at the size of 2.
   * \return The number of points found.
   */
  int points_at_x(const Curve_2& arc, const Point_2& p, Point_2* ps) const {
    // Get the y coordinates of the points on the conic.
    Algebraic ys[2];
    int n = conic_get_y_coordinates(arc, p.x(), ys);

    // Find all the points that are contained in the arc.
    int m = 0;

    for (int i = 0; i < n; ++i) {
      ps[m] = Point_2(p.x(), ys[i]);
      if (arc.is_full_conic() || is_between_endpoints(arc, ps[m])) ++m;
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
  int points_at_y(const Curve_2& arc, const Point_2& p, Point_2* ps) const {
    // Get the y coordinates of the points on the conic.
    Algebraic xs[2];
    int n = conic_get_x_coordinates(arc, p.y(), xs);

    // Find all the points that are contained in the arc.
    int m = 0;

    for (int i = 0; i < n; ++i) {
      ps[m] = Point_2(xs[i], p.y());
      if (arc.is_full_conic() || is_between_endpoints(arc, ps[m])) ++m;
    }

    // Return the number of points on the arc.
    CGAL_assertion(m <= 2);
    return m;
  }
};

#include <CGAL/enable_warnings.h>

} //namespace CGAL
#endif
