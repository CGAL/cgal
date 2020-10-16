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

#include <CGAL/atomic.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_geometry_traits/Conic_arc_2.h>
#include <CGAL/Arr_geometry_traits/Conic_x_monotone_arc_2.h>
#include <CGAL/Arr_geometry_traits/Conic_point_2.h>

namespace CGAL {

/*!
 * \class A traits class for maintaining an arrangement of conic arcs (bounded
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
template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class Arr_conic_traits_2
{
public:

  typedef Rat_kernel_                     Rat_kernel;
  typedef Alg_kernel_                     Alg_kernel;
  typedef Nt_traits_                      Nt_traits;

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

  mutable Intersection_map  inter_map;  // Mapping conic pairs to their
                                        // intersection points.

public:

  /*!
   * Default constructor.
   */
  Arr_conic_traits_2 ()
  {}

  /*! Get the next conic index. */
  static unsigned int get_index ()
  {
#ifdef CGAL_NO_ATOMIC
    static unsigned int index;
#else
    static CGAL::cpp11::atomic<unsigned int> index;
#endif
    return (++index);
  }

  /// \name Basic functor definitions.
  //@{

  class Compare_x_2
  {
  public:
    /*!
     * Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) < x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    Comparison_result operator() (const Point_2 & p1, const Point_2 & p2) const
    {
      Alg_kernel   ker;
      return (ker.compare_x_2_object() (p1, p2));
    }
  };

  /*! Get a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object () const
  {
    return Compare_x_2();
  }

  class Compare_xy_2
  {
  public:
    /*!
     * Compares two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2);
     *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL if the two points are equal.
     */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      Alg_kernel   ker;
      return (ker.compare_xy_2_object() (p1, p2));
    }
  };

  /*! Get a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return Compare_xy_2();
  }

  class Construct_min_vertex_2
  {
  public:
    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The left endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2 & cv) const
    {
      return (cv.left());
    }
  };

  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return Construct_min_vertex_2();
  }

  class Construct_max_vertex_2
  {
  public:
    /*!
     * Get the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \return The right endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2 & cv) const
    {
      return (cv.right());
    }
  };

  /*! Get a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2();
  }

  class Is_vertical_2
  {
  public:
    /*!
     * Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv) const
    {
      return (cv.is_vertical());
    }
  };

  /*! Get an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object () const
  {
    return Is_vertical_2();
  }

  class Compare_y_at_x_2
  {
  public:
    /*!
     * Return the location of the given point with respect to the input curve.
     * \param cv The curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator() (const Point_2 & p,
                                  const X_monotone_curve_2 & cv) const
    {
      Alg_kernel   ker;

      if (cv.is_vertical())
      {
        // A special treatment for vertical segments:
        // In case p has the same x c-ordinate of the vertical segment, compare
        // it to the segment endpoints to determine its position.
        Comparison_result res1 = ker.compare_y_2_object() (p, cv.left());
        Comparison_result res2 = ker.compare_y_2_object() (p, cv.right());

        if (res1 == res2)
          return (res1);
        else
          return (EQUAL);
      }

      // Check whether the point is exactly on the curve.
      if (cv.contains_point(p))
        return (EQUAL);

      // Get a point q on the x-monotone arc with the same x coordinate as p.
      Comparison_result  x_res;
      Point_2            q;

      if ((x_res = ker.compare_x_2_object() (p, cv.left())) == EQUAL)
      {
        q = cv.left();
      }
      else
      {
        CGAL_precondition (x_res != SMALLER);

        if ((x_res = ker.compare_x_2_object() (p, cv.right())) == EQUAL)
        {
          q = cv.right();
        }
        else
        {
          CGAL_precondition (x_res != LARGER);

          q = cv.point_at_x (p);
        }
      }

      // Compare p with the a point of the curve with the same x coordinate.
      return (ker.compare_y_2_object() (p, q));
    }
  };

  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return Compare_y_at_x_2();
  }

  class Compare_y_at_x_left_2
  {
  public:
    /*!
     * Compares the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const
    {
      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition (cv1.contains_point (p) &&
                         cv2.contains_point (p));

      CGAL_precondition_code (
        Alg_kernel   ker;
      );
      CGAL_precondition (ker.compare_xy_2_object() (p,
                                                    cv1.left()) == LARGER &&
                         ker.compare_xy_2_object() (p,
                                                    cv2.left()) == LARGER);

      // If one of the curves is vertical, it is below the other one.
      if (cv1.is_vertical())
      {
        if (cv2.is_vertical())
          // Both are vertical:
          return (EQUAL);
        else
          return (SMALLER);
      }
      else if (cv2.is_vertical())
      {
        return (LARGER);
      }

      // Compare the two curves immediately to the left of p:
      return (cv1.compare_to_left (cv2, p));
    }
  };

  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object () const
  {
    return Compare_y_at_x_left_2();
  }

  class Compare_y_at_x_right_2
  {
  public:
    /*!
     * Compares the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const
    {
      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition (cv1.contains_point (p) &&
                         cv2.contains_point (p));

      CGAL_precondition_code (
        Alg_kernel   ker;
      );

      CGAL_precondition (ker.compare_xy_2_object() (p,
                                                    cv1.right()) == SMALLER &&
                         ker.compare_xy_2_object() (p,
                                                    cv2.right()) == SMALLER);

      // If one of the curves is vertical, it is above the other one.
      if (cv1.is_vertical())
      {
        if (cv2.is_vertical())
          // Both are vertical:
          return (EQUAL);
        else
          return (LARGER);
      }
      else if (cv2.is_vertical())
      {
        return (SMALLER);
      }

      // Compare the two curves immediately to the right of p:
      return (cv1.compare_to_right (cv2, p));
    }
  };

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return Compare_y_at_x_right_2();
  }

  class Equal_2
  {
  public:
    /*!
     * Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      if (&cv1 == &cv2)
        return (true);

      return (cv1.equals (cv2));
    }

    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator() (const Point_2& p1, const Point_2& p2) const
    {
      if (&p1 == &p2)
        return (true);

      Alg_kernel   ker;
      return (ker.compare_xy_2_object() (p1, p2) == EQUAL);
    }
  };

  /*! Get an Equal_2 functor object. */
  Equal_2 equal_2_object () const
  {
    return Equal_2();
  }
  //@}

  /// \name Intersections, subdivisions, and mergings
  //@{

  /*! \class Make_x_monotone_2
   * A functor for subdividing curves into x-monotone curves.
   */
  class Make_x_monotone_2 {
    typedef Arr_conic_traits_2 <Rat_kernel_, Alg_kernel_, Nt_traits_>    Self;

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
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    {
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
          Alg_kernel_ ker;
          typename Alg_kernel_::Line_2 line =
            ker.construct_line_2_object()(vtan_ps[0], vtan_ps[1]);
          const Comparison_result start_pos =
            ker.compare_y_at_x_2_object() (cv.source(), line);
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

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(); }

  class Split_2
  {
  public:
    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint).
     * \param c2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator() (const X_monotone_curve_2& cv, const Point_2 & p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      cv.split (p, c1, c2);
      return;
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2();
  }

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
    { return (cv1.intersect(cv2, _inter_map, oi)); }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const { return (Intersect_2(inter_map)); }

  class Are_mergeable_2
  {
  public:
    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      return (cv1.can_merge_with (cv2));
    }
  };

  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const
  {
    return Are_mergeable_2();
  }

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2
  {
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
    /*!
     * Merge two given x-monotone curves into a single curve (segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable.
     */
    void operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c) const
    {
      CGAL_precondition(m_traits->are_mergeable_2_object()(cv2, cv1));

      c = cv1;
      c.merge (cv2);
    }
  };

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object() const
  {
    return Merge_2(this);
  }

  //@}

  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  typedef double                          Approximate_number_type;

  class Approximate_2
  {
  public:

    /*!
     * Return an approximation of a point coordinate.
     * \param p The exact point.
     * \param i The coordinate index (either 0 or 1).
     * \pre i is either 0 or 1.
     * \return An approximation of p's x-coordinate (if i == 0), or an
     *         approximation of p's y-coordinate (if i == 1).
     */
    Approximate_number_type operator() (const Point_2& p,
                                        int i) const
    {
      CGAL_precondition (i == 0 || i == 1);

      if (i == 0)
        return (CGAL::to_double(p.x()));
      else
        return (CGAL::to_double(p.y()));
    }
  };

  /*! Get an Approximate_2 functor object. */
  Approximate_2 approximate_2_object () const
  {
    return Approximate_2();
  }

  class Construct_x_monotone_curve_2
  {
  public:

    /*!
     * Return an x-monotone curve connecting the two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q must not be the same.
     * \return A segment connecting p and q.
     */
    X_monotone_curve_2 operator() (const Point_2& p,
                                   const Point_2& q) const
    {
      return (X_monotone_curve_2 (p, q));
    }
  };

  /*! Get a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object () const
  {
    return Construct_x_monotone_curve_2();
  }
  //@}

  /// \name Functor definitions for the Boolean set-operation traits.
  //@{

  class Compare_endpoints_xy_2
  {
  public:

    /*!
     * Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param cv The curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv) const
    {
      if (cv.is_directed_right())
        return (SMALLER);
      else
        return (LARGER);
    }
  };

  /*! Get a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  {
    return Compare_endpoints_xy_2();
  }

  class Construct_opposite_2
  {
  public:

    /*!
     * Construct an opposite x-monotone (with swapped source and target).
     * \param cv The curve.
     * \return The opposite curve.
     */
    X_monotone_curve_2 operator() (const X_monotone_curve_2& cv) const
    {
      return (cv.flip());
    }
  };

  /*! Get a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  {
    return Construct_opposite_2();
  }

  class Trim_2
  {
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
};

#include <CGAL/enable_warnings.h>

} //namespace CGAL
#endif
