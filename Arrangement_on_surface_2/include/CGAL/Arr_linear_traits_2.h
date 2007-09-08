// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_LINEAR_TRAITS_2_H
#define CGAL_ARR_LINEAR_TRAITS_2_H

/*! \file
 * The traits-class for handling linear objects (lines, rays and segments)
 * in the arrangement package.
 */

#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/intersections.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_geometry_traits/Segment_assertions.h>
#include <fstream>

CGAL_BEGIN_NAMESPACE

template <class Kernel_> class Arr_linear_object_2;

/*! \class
 * A traits class for maintaining an arrangement of linear objects (lines,
 * rays and segments), aoviding cascading of computations as much as possible.
 */
template <class Kernel_>
class Arr_linear_traits_2
{
  friend class Arr_linear_object_2<Kernel_>;

public:

  typedef Kernel_                         Kernel;
  typedef typename Kernel::FT             FT;

  typedef typename Algebraic_structure_traits<FT>::Is_exact 
                                          Has_exact_division;

  // Category tags:
  typedef Tag_true                        Has_left_category;
  typedef Tag_true                        Has_merge_category;
  typedef Tag_true                        Has_boundary_category;

  typedef typename Kernel::Line_2         Line_2;
  typedef typename Kernel::Ray_2          Ray_2;
  typedef typename Kernel::Segment_2      Segment_2;

  typedef CGAL::Segment_assertions<Arr_linear_traits_2<Kernel> >
                                          Segment_assertions;

  /*!
   * \class Representation of a linear with cached data.
   */
  class _Linear_object_cached_2
  {
  public:

    typedef typename Kernel::Line_2                Line_2;
    typedef typename Kernel::Ray_2                 Ray_2;
    typedef typename Kernel::Segment_2             Segment_2;
    typedef typename Kernel::Point_2               Point_2;

  protected:

    Line_2    l;                // The supporting line.
    Point_2   ps;               // The source point (if exists).
    Point_2   pt;               // The target point (if exists).
    bool      has_source;       // Is the source point valid
                                // (false for a line).
    bool      has_target;       // Is the target point valid
                                // (false for a line and for a ray).
    bool      is_right;         // Is the object directed to the right
                                // (for segments and rays).
    bool      is_vert;          // Is this a vertical object.
    bool      is_horiz;         // Is this a horizontal object.
    bool      has_pos_slope;    // Does the supporting line has a positive
                                // slope (if all three flags is_vert, is_horiz
                                // and has_pos_slope are false, then the line
                                // has a negative slope).
    bool      is_degen;         // Is the object degenerate (a single point).

  public:

    /*!
     * Default constructor.
     */
    _Linear_object_cached_2 () :
      has_source (true),
      has_target (true),
      is_vert (false),
      is_horiz (false),
      has_pos_slope (false),
      is_degen (true)
    {}

    /*!
     * Constructor for segment from two points.
     * \param p1 source point.
     * \param p2 target point.
     * \pre The two points must not be equal.
     */
    _Linear_object_cached_2(const Point_2& source, const Point_2& target) :
      ps (source),
      pt (target),
      has_source (true),
      has_target (true)
    {
      Kernel   kernel;

      Comparison_result  res = kernel.compare_xy_2_object()(source, target);
      is_degen = (res == EQUAL);
      is_right = (res == SMALLER);

      CGAL_precondition_msg (! is_degen,
                             "Cannot construct a degenerate segment.");

      l = kernel.construct_line_2_object()(source, target);
      is_vert = kernel.is_vertical_2_object()(l);
      is_horiz = kernel.is_horizontal_2_object()(l);
      has_pos_slope = _has_positive_slope();
    }

    /*!
     * Constructor from a segment.
     * \param seg The segment.
     * \pre The segment is not degenerate.
     */
    _Linear_object_cached_2 (const Segment_2& seg)
    {
      Kernel   kernel;

      CGAL_assertion_msg (! kernel.is_degenerate_2_object() (seg),
                          "Cannot construct a degenerate segment.");

      typename Kernel_::Construct_vertex_2
        construct_vertex = kernel.construct_vertex_2_object();

      ps = construct_vertex(seg, 0);
      has_source = true;
      pt = construct_vertex(seg, 1);
      has_target = true;

      Comparison_result  res = kernel.compare_xy_2_object()(ps, pt);
      
      CGAL_assertion (res != EQUAL);
      is_degen = false;
      is_right = (res == SMALLER);

      l = kernel.construct_line_2_object()(seg);
      is_vert = kernel.is_vertical_2_object()(seg);
      is_horiz = kernel.is_horizontal_2_object()(seg);
      has_pos_slope = _has_positive_slope();
    }

    /*!
     * Constructor from a ray.
     * \param ray The ray.
     * \pre The ray is not degenerate.
     */
    _Linear_object_cached_2 (const Ray_2& ray)
    {
      Kernel   kernel;

      CGAL_assertion_msg (! kernel.is_degenerate_2_object() (ray),
                          "Cannot construct a degenerate ray.");

      typename Kernel_::Construct_point_on_2
        construct_vertex = kernel.construct_point_on_2_object();

      ps = construct_vertex(ray, 0);         // The source point.
      has_source = true;
      pt = construct_vertex(ray, 1);         // Some point on the ray.
      has_target = false;

      Comparison_result  res = kernel.compare_xy_2_object()(ps, pt);
      CGAL_assertion (res != EQUAL);
      is_degen = false;
      is_right = (res == SMALLER);

      l = kernel.construct_line_2_object()(ray);
      is_vert = kernel.is_vertical_2_object()(ray);
      is_horiz = kernel.is_horizontal_2_object()(ray);
      has_pos_slope = _has_positive_slope();
    }

    /*!
     * Constructor from a line.
     * \param ln The line.
     * \pre The line is not degenerate.
     */
    _Linear_object_cached_2 (const Line_2& ln) :
      l (ln),
      has_source (false),
      has_target (false)
    {
      Kernel   kernel;

      CGAL_assertion_msg (! kernel.is_degenerate_2_object() (ln),
                          "Cannot construct a degenerate line.");

      typename Kernel_::Construct_point_on_2
        construct_vertex = kernel.construct_point_on_2_object();

      ps = construct_vertex(ln, 0);         // Some point on the line.
      has_source = false;
      pt = construct_vertex(ln, 1);         // Some point further on the line.
      has_target = false;

      Comparison_result  res = kernel.compare_xy_2_object()(ps, pt);
      CGAL_assertion (res != EQUAL);
      is_degen = false;
      is_right = (res == SMALLER);

      is_vert = kernel.is_vertical_2_object()(ln);
      is_horiz = kernel.is_horizontal_2_object()(ln);
      has_pos_slope = _has_positive_slope();
    }

    /*!
     * Check if the x-coordinate of the left point is infinite.
     * \return MINUS_INFINITY if the left point is at x = -oo;
     *         NO_BOUNDARY if the x-coordinate is finite.
     */
    Boundary_type left_infinite_in_x () const
    {
      if (is_vert || is_degen)
        return (NO_BOUNDARY);

      if (is_right)
        return (has_source ? NO_BOUNDARY : MINUS_INFINITY);
      else
        return (has_target ? NO_BOUNDARY : MINUS_INFINITY);
    }

    /*!
     * Check if the y-coordinate of the left point is infinite.
     * \return MINUS_INFINITY if the left point is at y = -oo;
     *         PLUS_INFINITY if the left point is at y = +oo;
     *         NO_BOUNDARY if the y-coordinate is finite.
     */
    Boundary_type left_infinite_in_y () const
    {
      if (is_horiz || is_degen)
        return (NO_BOUNDARY);

      if (is_vert)
      {
        if (is_right)
          return (has_source ? NO_BOUNDARY : MINUS_INFINITY);
        else
          return (has_target ? NO_BOUNDARY : MINUS_INFINITY);
      }

      if ((is_right && has_source) || (! is_right && has_target))
          return (NO_BOUNDARY);

      return (has_pos_slope ? MINUS_INFINITY : PLUS_INFINITY);
    }

    /*!
     * Check if the left point is finite.
     */
    bool has_left () const
    {
      if (is_right)
        return (has_source);
      else 
        return (has_target);
    }

    /*!
     * Get the (lexicographically) left endpoint.
     * \pre The left point is finite.
     */
    const Point_2& left () const
    {
      CGAL_precondition (has_left());
      return (is_right ? ps : pt);
    }

    /*!
     * Set the (lexicographically) left endpoint.
     * \param p The point to set.
     * \pre p lies on the supporting line to the left of the right endpoint.
     */
    void set_left (const Point_2& p, bool check_validity = true)
    {
      CGAL_precondition (! is_degen);

      CGAL_precondition_code (
        Kernel    kernel;
      );
      CGAL_precondition 
        (Segment_assertions::_assert_is_point_on (p, l, 
                                                  Has_exact_division()) &&
         (! check_validity || ! has_right() ||
          kernel.compare_xy_2_object() (p, right()) == SMALLER));

      if (is_right)
      {
        ps = p;
        has_source = true;
      }
      else
      {
        pt = p;
        has_target = true;
      }
    }

    /*!
     * Set the (lexicographically) left endpoint as infinite.
     */
    void set_left ()
    {
      CGAL_precondition (! is_degen);

      if (is_right)
        has_source = false;
      else
        has_target = false;
    }

    /*!
     * Check if the x-coordinate of the right point is infinite.
     * \return PLUS_INFINITY if the left point is at x = +oo;
     *         NO_BOUNDARY if the x-coordinate is finite.
     */
    Boundary_type right_infinite_in_x () const
    {
      if (is_vert || is_degen)
        return (NO_BOUNDARY);

      if (is_right)
        return (has_target ? NO_BOUNDARY : PLUS_INFINITY);
      else
        return (has_source ? NO_BOUNDARY : PLUS_INFINITY);
    }

    /*!
     * Check if the y-coordinate of the right point is infinite.
     * \return MINUS_INFINITY if the right point is at y = -oo;
     *         PLUS_INFINITY if the right point is at y = +oo;
     *         NO_BOUNDARY if the y-coordinate is finite.
     */
    Boundary_type right_infinite_in_y () const
    {
      if (is_horiz || is_degen)
        return (NO_BOUNDARY);

      if (is_vert)
      {
        if (is_right)
          return (has_target ? NO_BOUNDARY : PLUS_INFINITY);
        else
          return (has_source ? NO_BOUNDARY : PLUS_INFINITY);
      }

      if ((is_right && has_target) || (! is_right && has_source))
          return (NO_BOUNDARY);

      return (has_pos_slope ? PLUS_INFINITY : MINUS_INFINITY);
    }

    /*!
     * Check if the right point is finite.
     */
    bool has_right () const
    {
      if (is_right)
        return (has_target);
      else 
        return (has_source);
    }

    /*!
     * Get the (lexicographically) right endpoint.
     * \pre The right endpoint is finite.
     */
    const Point_2& right () const
    {
      CGAL_precondition (has_right());
      return (is_right ? pt : ps);
    }

    /*!
     * Set the (lexicographically) right endpoint.
     * \param p The point to set.
     * \pre p lies on the supporting line to the right of the left endpoint.
     */
    void set_right (const Point_2& p, bool check_validity = true)
    {
      CGAL_precondition (! is_degen);
      CGAL_precondition_code (
        Kernel    kernel;
      );
      CGAL_precondition 
        (Segment_assertions::_assert_is_point_on (p, l, 
                                                  Has_exact_division()) &&
         (! check_validity || ! has_left() ||
          kernel.compare_xy_2_object() (p, left()) == LARGER));

      if (is_right)
      {
        pt = p;
        has_target = true;
      }
      else
      {
        ps = p;
        has_source = true;
      }
    }

    /*!
     * Set the (lexicographically) right endpoint as infinite.
     */
    void set_right ()
    {
      CGAL_precondition (! is_degen);

      if (is_right)
        has_target = false;
      else
        has_source = false;
    }

    /*!
     * Get the supporting line.
     */
    const Line_2& supp_line () const
    {
      CGAL_precondition (! is_degen);
      return (l);
    }

    /*!
     * Check if the curve is vertical.
     */
    bool is_vertical () const
    {
      CGAL_precondition (! is_degen);
      return (is_vert);
    }

    /*!
     * Check if the curve is degenerate.
     */
    bool is_degenerate () const
    {
      return (is_degen);
    }

    /*!
     * Check if the curve is directed lexicographic from left to right
     */
    bool is_directed_right () const
    {
      return (is_right);
    }

    /*!
     * Check if the given point is in the x-range of the object.
     * \param p The query point.
     * \return (true) is in the x-range of the segment; (false) if it is not.
     */
    bool is_in_x_range (const Point_2& p) const
    {
      Kernel                          kernel;
      typename Kernel_::Compare_x_2   compare_x = kernel.compare_x_2_object();
      Comparison_result               res1;

      if (left_infinite_in_x() == NO_BOUNDARY)
      {
        if (left_infinite_in_y() != NO_BOUNDARY)
          // Compare with some point on the curve.
          res1 = compare_x (p, ps);
        else
          res1 = compare_x (p, left());
      }
      else
      {
        // p is obviously to the right.
        res1 = LARGER;
      }

      if (res1 == SMALLER)
        return (false);
      else if (res1 == EQUAL)
        return (true);

      Comparison_result               res2;

      if (right_infinite_in_x() == NO_BOUNDARY)
      {
        if (right_infinite_in_y() != NO_BOUNDARY)
          // Compare with some point on the curve.
          res2 = compare_x (p, ps);
        else
          res2 = compare_x (p, right());
      }
      else
      {
        // p is obviously to the right.
        res2 = SMALLER;          
      }

      return (res2 != LARGER);
    }

    /*!
     * Check if the given point is in the y-range of the object.
     * \param p The query point.
     * \pre The object is vertical.
     * \return (true) is in the y-range of the segment; (false) if it is not.
     */
    bool is_in_y_range (const Point_2& p) const
    {
      CGAL_precondition (is_vertical());

      Kernel                          kernel;
      typename Kernel_::Compare_y_2   compare_y = kernel.compare_y_2_object();
      Boundary_type                   inf = left_infinite_in_y();
      Comparison_result               res1;

      CGAL_assertion (inf != PLUS_INFINITY);
      if (inf == NO_BOUNDARY)
        res1 = compare_y (p, left());
      else
        res1 = LARGER;           // p is obviously above.

      if (res1 == SMALLER)
        return (false);
      else if (res1 == EQUAL)
        return (true);

      Comparison_result               res2;

      inf = right_infinite_in_y();
      CGAL_assertion (inf != MINUS_INFINITY);
      if (inf == NO_BOUNDARY)
        res2 = compare_y (p, right());
      else
        res2 = SMALLER;          // p is obviously below.

      return (res2 != LARGER);
    }

  private:

    /*!
     * Determine if the supporting line has a positive slope.
     */
    bool _has_positive_slope () const
    {
      if (is_vert)
        return (true);

      if (is_horiz)
        return (false);

      // Construct a horizontal line and compare its slope the that of l.
      Kernel     kernel;
      Line_2     l_horiz = kernel.construct_line_2_object() (Point_2 (0, 0),
                                                             Point_2 (1, 0));

      return (kernel.compare_slope_2_object() (l, l_horiz) == LARGER);
    }
  };

public:

  // Traits objects
  typedef typename Kernel::Point_2              Point_2;
  typedef Arr_linear_object_2<Kernel>           X_monotone_curve_2;
  typedef Arr_linear_object_2<Kernel>           Curve_2;

public:

  /*!
   * Default constructor.
   */
  Arr_linear_traits_2 ()
  {}

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
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      Kernel    kernel;

      return (kernel.compare_x_2_object()(p1, p2));
    }

    /*!
     * Compare the relative positions of a vertical curve and another given
     * curves at y = +/- oo.
     * \param p A reference point; we refer to a vertical line incident to p.
     * \param cv The compared curve.
     * \param ind MIN_END if we refer to cv's minimal end,
     *            MIN_END if we refer to its maximal end.
     * \pre cv's relevant end is defined at y = +/- oo.
     * \return SMALLER if p lies to the left of cv;
     *         LARGER if p lies to the right cv;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv,
                                  Curve_end) const
    {
      CGAL_precondition (! cv.is_degenerate());
      CGAL_precondition (cv.is_vertical());

      Kernel                    kernel;
      return (kernel.compare_x_at_y_2_object() (p, cv.supp_line()));
    }

    /*!
     * Compare the relative positions of two curves at y = +/- oo.
     * \param cv1 The first curve.
     * \param ind1 MIN_END if we refer to cv1's minimal end,
     *             MIN_END if we refer to its maximal end.
     * \param cv2 The second curve.
     * \param ind2 MIN_END if we refer to cv2's minimal end,
     *             MIN_END if we refer to its maximal end.
     * \pre The curves are defined at y = +/- oo.
     * \return SMALLER if cv1 lies to the left of cv2;
     *         LARGER if cv1 lies to the right cv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  Curve_end /* ind1 */,
                                  const X_monotone_curve_2& cv2,
                                  Curve_end /* ind2 */) const
    {
      CGAL_precondition (! cv1.is_degenerate());
      CGAL_precondition (! cv2.is_degenerate());
      CGAL_precondition (cv1.is_vertical());
      CGAL_precondition (cv2.is_vertical());

      Kernel                    kernel;
      typename Kernel::Point_2 p = kernel.construct_point_2_object() (ORIGIN);
      return (kernel.compare_x_at_y_2_object() (p,
                                                cv1.supp_line(),
                                                cv2.supp_line()));
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
     * Compare two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2);
     *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL if the two points are equal.
     */
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      Kernel    kernel;
      return (kernel.compare_xy_2_object()(p1, p2));
    }
  };

  /*! Get a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return Compare_xy_2();
  }

  class Boundary_in_x_2
  {
  public:
    /*!
     * Check if an end of a given x-monotone curve is infinite at x.
     * \param cv The curve.
     * \param ind MIN_END if we refer to cv's minimal end,
     *            MAX_END if we refer to its maximal end.
     * \return MINUS_INFINITY if the curve end lies at x = -oo;
     *         NO_BOUNDARY if the curve end has a finite x-coordinate;
     *         PLUS_INFINITY if the curve end lies at x = +oo.
     */
    Boundary_type operator() (const X_monotone_curve_2& cv,
                              Curve_end ind) const
    {
      CGAL_precondition (! cv.is_degenerate());

      if (ind == MIN_END)
        return (cv.left_infinite_in_x());
      else
        return (cv.right_infinite_in_x());
    }
  };

  /*! Get an Boundary_in_x_2 functor object. */
  Boundary_in_x_2 boundary_in_x_2_object () const
  {
    return Boundary_in_x_2();
  }

  class Boundary_in_y_2
  {
  public:
    /*!
     * Check if an end of a given x-monotone curve is infinite at y.
     * \param cv The curve.
     * \param ind MIN_END if we refer to cv's minimal end,
     *            MAX_END if we refer to its maximal end.
     * \return MINUS_INFINITY if the curve end lies at y = -oo;
     *         NO_BOUNDARY if the curve end has a finite y-coordinate;
     *         PLUS_INFINITY if the curve end lies at y = +oo.
     */
    Boundary_type operator() (const X_monotone_curve_2& cv,
                              Curve_end ind) const
    {
      CGAL_precondition (! cv.is_degenerate());

      if (ind == MIN_END)
        return (cv.left_infinite_in_y());
      else
        return (cv.right_infinite_in_y());
    }
  };

  /*! Get an Boundary_in_y_2 functor object. */
  Boundary_in_y_2 boundary_in_y_2_object () const
  {
    return Boundary_in_y_2();
  }

  class Construct_min_vertex_2
  {
  public:
    /*!
     * Get the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \pre The left end of cv is a valid (bounded) point.
     * \return The left endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2& cv) const
    {
      CGAL_precondition (! cv.is_degenerate());
      CGAL_precondition (cv.has_left());

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
     * \pre The right end of cv is a valid (bounded) point.
     * \return The right endpoint.
     */
    const Point_2& operator() (const X_monotone_curve_2& cv) const
    {
      CGAL_precondition (! cv.is_degenerate());
      CGAL_precondition (cv.has_right());

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
      CGAL_precondition (! cv.is_degenerate());
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
    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv) const
    {
      CGAL_precondition (! cv.is_degenerate());
      CGAL_precondition (cv.is_in_x_range (p));

      Kernel    kernel;

      if (! cv.is_vertical())
      {
        // Compare p with the segment's supporting line.
        return (kernel.compare_y_at_x_2_object()(p, cv.supp_line()));
      }
      else
      {
        // Compare with the vertical segment's end-points.
        typename Kernel::Compare_y_2  compare_y = kernel.compare_y_2_object();
        const Comparison_result res1 =
          cv.has_left() ? compare_y (p, cv.left()) : LARGER;
        const Comparison_result res2 = 
          cv.has_right() ? compare_y (p, cv.right()) : SMALLER;

        if (res1 == res2)
          return (res1);
        else
          return (EQUAL);
      }
    }

    /*!
     * Compare the relative y-positions of two curves at x = +/- oo.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param ind MIN_END if we compare at x = -oo;
     *            MAX_END if we compare at x = +oo.
     * \pre The curves are defined at x = +/- oo.
     * \return SMALLER if cv1 lies below cv2;
     *         LARGER if cv1 lies above cv2;
     *         EQUAL in case of an overlap.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2, 
                                  Curve_end ind) const
    {
      // Make sure both curves are defined at x = -oo (or at x = +oo).
      CGAL_precondition (! cv1.is_degenerate());
      CGAL_precondition (! cv2.is_degenerate());
      CGAL_precondition ((ind == MIN_END &&
                          cv1.left_infinite_in_x() == MINUS_INFINITY &&
                          cv2.left_infinite_in_x() == MINUS_INFINITY) ||
                         (ind == MAX_END &&
                          cv1.right_infinite_in_x() == PLUS_INFINITY &&
                          cv2.right_infinite_in_x() == PLUS_INFINITY));

      // Compare the slopes of the two supporting lines.
      Kernel                    kernel;
      const Comparison_result   res_slopes =
        kernel.compare_slope_2_object() (cv1.supp_line(),
                                         cv2.supp_line());

      if (res_slopes == EQUAL)
      {
        // In case the two supporting line are parallel, compare their
        // relative position at x = 0, which is the same as their position
        // at infinity.
        typename Kernel::Point_2 p = kernel.construct_point_2_object() (ORIGIN);
        return (kernel.compare_y_at_x_2_object() (p,
                                                  cv1.supp_line(),
                                                  cv2.supp_line()));
      }

      if (ind == MIN_END)
        // Flip the slope result if we compare at x = -oo:
        return ((res_slopes == LARGER) ? SMALLER : LARGER);

      // If we compare at x = +oo, the slope result is what we need:
      return (res_slopes);
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
     * Compare the y value of two x-monotone curves immediately to the left
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
                                  const Point_2& CGAL_precondition_code(p) ) const
    {
      CGAL_precondition (! cv1.is_degenerate());
      CGAL_precondition (! cv2.is_degenerate());

      Kernel                        kernel;

      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition_code (
        typename Kernel::Compare_xy_2 compare_xy = 
                                                  kernel.compare_xy_2_object();
      );

      CGAL_precondition 
        (Segment_assertions::_assert_is_point_on (p, cv1, 
                                                  Has_exact_division()) &&
         Segment_assertions::_assert_is_point_on (p, cv2,
                                                  Has_exact_division()));

      CGAL_precondition ((! cv1.has_left() ||
                          compare_xy(cv1.left(), p) == SMALLER) &&
                         (! cv2.has_left() ||
                          compare_xy(cv2.left(), p) == SMALLER));

      // Compare the slopes of the two segments to determine thir relative
      // position immediately to the left of q.
      // Notice we use the supporting lines in order to compare the slopes,
      // and that we swap the order of the curves in order to obtain the
      // correct result to the left of p.
      return (kernel.compare_slope_2_object()(cv2.supp_line(),
                                              cv1.supp_line()));
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
     * Compare the y value of two x-monotone curves immediately to the right
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
                                  const Point_2& CGAL_precondition_code(p) ) const
    {
      CGAL_precondition (! cv1.is_degenerate());
      CGAL_precondition (! cv2.is_degenerate());

      Kernel                        kernel;

      // Make sure that p lies on both curves, and that both are defined to its
      // right (so their right endpoint is lexicographically larger than p).
      CGAL_precondition_code (
        typename Kernel::Compare_xy_2 compare_xy = 
                                                 kernel.compare_xy_2_object();
      );

      CGAL_precondition
        (Segment_assertions::_assert_is_point_on (p, cv1, 
                                                  Has_exact_division()) &&
         Segment_assertions::_assert_is_point_on (p, cv2,
                                                  Has_exact_division()));

      CGAL_precondition ((! cv1.has_right() ||
                          compare_xy(cv1.right(), p) == LARGER) &&
                         (! cv2.has_right() ||
                          compare_xy(cv2.right(), p) == LARGER));

      // Compare the slopes of the two segments to determine thir relative
      // position immediately to the left of q.
      // Notice we use the supporting lines in order to compare the slopes.
      return (kernel.compare_slope_2_object()(cv1.supp_line(),
                                              cv2.supp_line()));
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
      CGAL_precondition (! cv1.is_degenerate());
      CGAL_precondition (! cv2.is_degenerate());

      Kernel                    kernel;
      typename Kernel::Equal_2  equal = kernel.equal_2_object();

      // Check that the two supporting lines are the same.
      if (! equal (cv1.supp_line(),
                   cv2.supp_line()) &&
          ! equal (cv1.supp_line(), 
                   kernel.construct_opposite_line_2_object()(cv2.supp_line())))
      {
        return (false);
      }

      // Check that either the two left endpoints are at infinity, or they
      // are bounded and equal.
      if ((cv1.has_left() != cv2.has_left()) ||
          (cv1.has_left() && ! equal (cv1.left(), cv2.left())))
      {
        return (false);
      }

      // Check that either the two right endpoints are at infinity, or they
      // are bounded and equal.
      return ((cv1.has_right() == cv2.has_right()) &&
              (! cv1.has_right() || equal (cv1.right(), cv2.right())));
    }

    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator() (const Point_2& p1, const Point_2& p2) const
    {
      Kernel    kernel;
      return (kernel.equal_2_object()(p1, p2));
    }
  };

  /*! Get an Equal_2 functor object. */
  Equal_2 equal_2_object () const
  {
    return Equal_2();
  }
  //@}

  /// \name Functor definitions for supporting intersections.
  //@{

  class Make_x_monotone_2
  {
  public:
    /*!
     * Cut the given curve into x-monotone subcurves and insert them into the
     * given output iterator. As segments are always x_monotone, only one
     * object will be contained in the iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The output
     *           object is a wrapper of an X_monotone_curve_2 which is
     *           essentially the same as the input curve.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi) const
    {
      // Wrap the curve with an object.
      *oi = make_object (cv);
      ++oi;

      return (oi);
    }
  };

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () const
  {
    return Make_x_monotone_2();
  }

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
    void operator() (const X_monotone_curve_2& cv, const Point_2& p,
                     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      CGAL_precondition (! cv.is_degenerate());

      // Make sure that p lies on the interior of the curve.
      CGAL_precondition_code (
        Kernel                        kernel;
        typename Kernel::Compare_xy_2 compare_xy = 
                                                 kernel.compare_xy_2_object();
      );

      CGAL_precondition
        (Segment_assertions::_assert_is_point_on (p, cv,
                                                  Has_exact_division()) &&
         (! cv.has_left() || compare_xy(cv.left(), p) == SMALLER) &&
         (! cv.has_right() || compare_xy(cv.right(), p) == LARGER));

      // Perform the split.
      c1 = cv;
      c1.set_right (p);

      c2 = cv;
      c2.set_left (p);

      return;
    }
  };

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2();
  }

  class Intersect_2
  {
  public:
    /*!
     * Find the intersections of the two given curves and insert them into the
     * given output iterator. As two segments may itersect only once, only a
     * single intersection will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi) const
    {
      CGAL_precondition (! cv1.is_degenerate());
      CGAL_precondition (! cv2.is_degenerate());

      // Intersect the two supporting lines.
      Kernel       kernel;
      CGAL::Object obj = kernel.intersect_2_object()(cv1.supp_line(),
                                                     cv2.supp_line());

      if (obj.is_empty())
      {
        // The supporting line are parallel lines and do not intersect:
        return (oi);
      }

      // Check if we have a single intersection point.
      const Point_2  *ip = object_cast<Point_2> (&obj);
      
      if (ip != NULL)
      {
        // Check if the intersection point ip lies on both segments.
        const bool    ip_on_cv1 = cv1.is_vertical() ? cv1.is_in_y_range(*ip) :
                                                      cv1.is_in_x_range(*ip);

        if (ip_on_cv1)
        {
          const bool  ip_on_cv2 = cv2.is_vertical() ? cv2.is_in_y_range(*ip) :
                                                      cv2.is_in_x_range(*ip);

          if (ip_on_cv2)
          {
            // Create a pair representing the point with its multiplicity,
            // which is always 1 for line segments.
            std::pair<Point_2, unsigned int>   ip_mult (*ip, 1);
            *oi = make_object (ip_mult);
            oi++;
          }
        }
        return (oi);
      }

      // In this case, the two supporting lines overlap.
      // We start with the entire cv1 curve as the overlapping subcurve,
      // then clip it to form the true overlapping curve.
      typename Kernel::Compare_xy_2  compare_xy = kernel.compare_xy_2_object();
      X_monotone_curve_2             ovlp = cv1;

      if (cv2.has_left())
      {
        // If the left endpoint of cv2 is to the right of cv1's left endpoint,
        // clip the overlapping subcurve.
        if (! cv1.has_left())
        {
          ovlp.set_left (cv2.left(), false);
        }
        else
        {
          if (compare_xy (cv1.left(), cv2.left()) == SMALLER)
            ovlp.set_left (cv2.left(), false);
        }
      }

      if (cv2.has_right())
      {
        // If the right endpoint of cv2 is to the left of cv1's right endpoint,
        // clip the overlapping subcurve.
        if (! cv1.has_right())
        {
          ovlp.set_right (cv2.right(), false);
        }
        else
        {
          if (compare_xy (cv1.right(), cv2.right()) == LARGER)
            ovlp.set_right (cv2.right(), false);
        }
      }

      // Examine the resulting subcurve.
      Comparison_result        res = SMALLER;

      if (ovlp.has_left() && ovlp.has_right())
        res = compare_xy (ovlp.left(), ovlp.right());

      if (res == SMALLER)
      {
        // We have discovered a true overlapping subcurve:
        *oi = make_object (ovlp);
        oi++;
      }
      else if (res == EQUAL)
      {
        // The two objects have the same supporting line, but they just share
        // a common endpoint. Thus we have an intersection point, but we leave
        // the multiplicity of this point undefined.
        std::pair<Point_2, unsigned int>   ip_mult (ovlp.left(), 0);
        *oi = make_object (ip_mult);
        oi++;
      }

      return (oi);
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return Intersect_2();
  }

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
      CGAL_precondition (! cv1.is_degenerate());
      CGAL_precondition (! cv2.is_degenerate());

      Kernel                    kernel;
      typename Kernel::Equal_2  equal = kernel.equal_2_object();

      // Check if the two curves have the same supporting line.
      if (! equal (cv1.supp_line(), cv2.supp_line()) && 
          ! equal (cv1.supp_line(), 
                   kernel.construct_opposite_line_2_object()(cv2.supp_line())))
        return (false);

      // Check if the left endpoint of one curve is the right endpoint of the
      // other.
      return ((cv1.has_right() && cv2.has_left() &&
               equal (cv1.right(), cv2.left())) ||
              (cv2.has_right() && cv1.has_left() &&
               equal (cv2.right(), cv1.left())));
    }
  };

  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const
  {
    return Are_mergeable_2();
  }

  class Merge_2
  {
  public:
    /*!
     * Merge two given x-monotone curves into a single curve (segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      same line and share a common endpoint.
     */
    void operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c) const
    {
      CGAL_precondition (! cv1.is_degenerate());
      CGAL_precondition (! cv2.is_degenerate());

      Kernel                    kernel;
      typename Kernel::Equal_2  equal = kernel.equal_2_object();

      CGAL_precondition
        (equal (cv1.supp_line(), 
                cv2.supp_line()) ||
         equal (cv1.supp_line(),
                kernel.construct_opposite_line_2_object()(cv2.supp_line())));

      // Check which curve extends to the right of the other.
      if (cv1.has_right() && cv2.has_left() &&
          equal (cv1.right(), cv2.left()))
      {
        // cv2 extends cv1 to the right.
        c = cv1;

        if (cv2.has_right())
          c.set_right (cv2.right());
        else
          c.set_right();      // Unbounded endpoint. 
      }
      else
      {
        CGAL_precondition (cv2.has_right() && cv1.has_left() &&
                           equal (cv2.right(), cv1.left()));

        // cv1 extends cv2 to the right.
        c = cv2;

        if (cv1.has_right())
          c.set_right (cv1.right());
        else
          c.set_right();      // Unbounded endpoint.
      }

      return;
    }
  };

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object () const
  {
    return Merge_2();
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
      Kernel     kernel;
      Segment_2  seg = kernel.construct_segment_2_object() (p, q);

      return (X_monotone_curve_2 (seg));
    }
  };

  /*! Get a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object () const
  {
    return Construct_x_monotone_curve_2();
  }
  //@}

};

/*!
 * \class A representation of a segment, as used by the Arr_segment_traits_2
 * traits-class.
 */
template <class Kernel_>
class Arr_linear_object_2 :
    public Arr_linear_traits_2<Kernel_>::_Linear_object_cached_2
{
  typedef typename Arr_linear_traits_2<Kernel_>::_Linear_object_cached_2 
                                                                   Base;

public:

  typedef Kernel_                                           Kernel;

  typedef typename Kernel::Point_2                          Point_2;
  typedef typename Kernel::Segment_2                        Segment_2;
  typedef typename Kernel::Ray_2                            Ray_2;
  typedef typename Kernel::Line_2                           Line_2;

public:

  /*!
   * Default constructor.
   */
  Arr_linear_object_2 () :
    Base()
  {}
    
  /*!
   * Constructor from two points.
   * \param s The source point.
   * \param t The target point.
   * \pre The two points must not be the same.
   */
  Arr_linear_object_2(const Point_2& s, const Point_2& t):
    Base(s, t)
  {}

  /*!
   * Constructor from a segment.
   * \param seg The segment.
   * \pre The segment is not degenerate.
   */
  Arr_linear_object_2 (const Segment_2& seg) :
    Base (seg)
  {}

  /*!
   * Constructor from a ray.
   * \param ray The segment.
   * \pre The ray is not degenerate.
   */
  Arr_linear_object_2 (const Ray_2& ray) :
    Base (ray)
  {}

  /*!
   * Constructor from a line.
   * \param line The line.
   * \pre The line is not degenerate.
   */
  Arr_linear_object_2 (const Line_2& line) :
    Base (line)
  {}

  /*!
   * Check if the object is actually a segment.
   */
  bool is_segment () const
  {
    return (! this->is_degen && this->has_source && this->has_target);
  }

  /*!
   * Cast to a segment.
   * \pre The linear object is really a segment.
   */
  Segment_2 segment () const
  {
    CGAL_precondition (is_segment());

    Kernel     kernel;
    Segment_2  seg = kernel.construct_segment_2_object() (this->ps, this->pt);
    return (seg);
  }

  /*!
   * Check if the object is actually a ray.
   */
  bool is_ray () const
  {
    return (! this->is_degen && 
            (this->has_source != this->has_target));
  }

  /*!
   * Cast to a ray.
   * \pre The linear object is really a ray.
   */
  Ray_2 ray () const
  {
    CGAL_precondition (is_ray());

    Kernel     kernel;
    Ray_2      ray;

    if (this->has_source)
      ray = kernel.construct_ray_2_object() (this->ps, this->l);
    else
      ray = kernel.construct_ray_2_object()
        (this->pt, 
         kernel.construct_opposite_line_2_object()(this->l));

    return (ray);
  }

  /*!
   * Check if the object is actually a line.
   */
  bool is_line () const
  {
    return (! this->is_degen && ! this->has_source && ! this->has_target);
  }

  /*!
   * Cast to a line.
   * \pre The linear object is really a line.
   */
  Line_2 line () const
  {
    CGAL_precondition (is_line());
    return (this->l);
  }

  /*!
   * Get the supporting line.
   * \pre The object is not a point.
   */
  const Line_2& supporting_line () const
  {
    CGAL_precondition (! this->is_degen);
    return (this->l);
  }

  /*!
   * Get the source point.
   * \pre The object is a point, a segment or a ray.
   */
  const Point_2& source() const
  {
    CGAL_precondition (! is_line());

    if (this->is_degen)
      return (this->ps);      // For a point.

    if (this->has_source)
      return (this->ps);      // For a segment or a ray.
    else
      return (this->pt);      // For a "flipped" ray.
  }

  /*!
   * Get the target point.
   * \pre The object is a point or a segment.
   */
  const Point_2& target() const
  {
    CGAL_precondition (! is_line() && ! is_ray());

    return (this->pt);
  }

  /*!
   * Create a bounding box for the linear object.
   */
  Bbox_2 bbox() const
  {
    CGAL_precondition(this->is_segment());
    Kernel     kernel;
    Segment_2  seg = kernel.construct_segment_2_object() (this->ps, this->pt);
    return (kernel.construct_bbox_2_object() (seg));
  }
};

/*!
 * Exporter for the segment class used by the traits-class.
 */
template <class Kernel, class OutputStream>
OutputStream& operator<< (OutputStream& os,
                          const Arr_linear_object_2<Kernel>& lobj)
{
  // Print a letter identifying the object type, then the object itself.
  if (lobj.is_segment())
    os << " S " << lobj.segment();
  else if (lobj.is_ray())
    os << " R " << lobj.ray();
  else
    os << " L " << lobj.line();

  return (os);
}

/*!
 * Importer for the segment class used by the traits-class.
 */
template <class Kernel, class InputStream>
InputStream& operator>> (InputStream& is, Arr_linear_object_2<Kernel>& lobj)
{
  // Read the object type.
  char        c;

  do
  {
    is >> c;
  } while ((c != 'S' && c != 's') &&
           (c != 'R' && c != 'r') &&
           (c != 'L' && c != 'l'));

  // Read the object accordingly.
  if (c == 'S' || c == 's')
  {
    typename Kernel::Segment_2  seg;
    is >> seg;
    lobj = seg;
  }
  else if (c == 'R' || c == 'r')
  {
    typename Kernel::Ray_2      ray;
    is >> ray;
    lobj = ray;
  }
  else
  {
    typename Kernel::Line_2     line;
    is >> line;
    lobj = line;
  }

  return (is);
}

CGAL_END_NAMESPACE

#endif
