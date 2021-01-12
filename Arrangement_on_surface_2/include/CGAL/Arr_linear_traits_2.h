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
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Waqar Khan        <wkhan@mpi-inf.mpg.de>
//            Efi fogel         <efifogel@gmail.com>

#ifndef CGAL_ARR_LINEAR_TRAITS_2_H
#define CGAL_ARR_LINEAR_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The traits-class for handling linear objects (lines, rays and segments)
 * in the arrangement package.
 */

#include <fstream>

#include <boost/variant.hpp>

#include <CGAL/tags.h>
#include <CGAL/intersections.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_geometry_traits/Segment_assertions.h>

namespace CGAL {

template <typename Kernel_> class Arr_linear_object_2;

/*! \class
 * A traits class for maintaining an arrangement of linear objects (lines,
 * rays and segments), aoviding cascading of computations as much as possible.
 */
template <typename Kernel_>
class Arr_linear_traits_2 : public Kernel_ {
  friend class Arr_linear_object_2<Kernel_>;

public:
  typedef Kernel_                         Kernel;
  typedef typename Kernel::FT             FT;

  typedef typename Algebraic_structure_traits<FT>::Is_exact
                                          Has_exact_division;

  // Category tags:
  typedef Tag_true                        Has_left_category;
  typedef Tag_true                        Has_merge_category;
  typedef Tag_false                       Has_do_intersect_category;

  typedef Arr_open_side_tag               Left_side_category;
  typedef Arr_open_side_tag               Bottom_side_category;
  typedef Arr_open_side_tag               Top_side_category;
  typedef Arr_open_side_tag               Right_side_category;

  typedef typename Kernel::Line_2         Line_2;
  typedef typename Kernel::Ray_2          Ray_2;
  typedef typename Kernel::Segment_2      Segment_2;

  typedef CGAL::Segment_assertions<Arr_linear_traits_2<Kernel> >
                                          Segment_assertions;

  /*!
   * \class Representation of a linear with cached data.
   */
  class _Linear_object_cached_2 {
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
    /*! Default constructor.
     */
    _Linear_object_cached_2() :
      has_source(true),
      has_target(true),
      is_vert(false),
      is_horiz(false),
      has_pos_slope(false),
      is_degen(true)
    {}

    /*! Constructor for segment from two points.
     * \param p1 source point.
     * \param p2 target point.
     * \pre The two points must not be equal.
     */
    _Linear_object_cached_2(const Point_2& source, const Point_2& target) :
      ps(source),
      pt(target),
      has_source(true),
      has_target(true)
    {
      Kernel kernel;

      Comparison_result res = kernel.compare_xy_2_object()(source, target);
      is_degen = (res == EQUAL);
      is_right = (res == SMALLER);

      CGAL_precondition_msg(! is_degen,
                            "Cannot construct a degenerate segment.");

      l = kernel.construct_line_2_object()(source, target);
      is_vert = kernel.is_vertical_2_object()(l);
      is_horiz = kernel.is_horizontal_2_object()(l);
      has_pos_slope = _has_positive_slope();
    }

    /*! Constructor from a segment.
     * \param seg The segment.
     * \pre The segment is not degenerate.
     */
    _Linear_object_cached_2(const Segment_2& seg) :
      has_source(true),
      has_target(true)
    {
      Kernel kernel;

      CGAL_assertion_msg(! kernel.is_degenerate_2_object()(seg),
                         "Cannot construct a degenerate segment.");

      auto construct_vertex = kernel.construct_vertex_2_object();
      ps = construct_vertex(seg, 0);
      pt = construct_vertex(seg, 1);

      Comparison_result res = kernel.compare_xy_2_object()(ps, pt);
      CGAL_assertion(res != EQUAL);
      is_degen = false;
      is_right = (res == SMALLER);

      l = kernel.construct_line_2_object()(seg);
      is_vert = kernel.is_vertical_2_object()(seg);
      is_horiz = kernel.is_horizontal_2_object()(seg);
      has_pos_slope = _has_positive_slope();
    }

    /*! Constructor from a ray.
     * \param ray The ray.
     * \pre The ray is not degenerate.
     */
    _Linear_object_cached_2(const Ray_2& ray) :
      has_source(true),
      has_target(false)
    {
      Kernel kernel;

      CGAL_assertion_msg(! kernel.is_degenerate_2_object()(ray),
                         "Cannot construct a degenerate ray.");

      auto construct_vertex = kernel.construct_point_on_2_object();
      ps = construct_vertex(ray, 0);         // The source point.
      pt = construct_vertex(ray, 1);         // Some point on the ray.

      Comparison_result  res = kernel.compare_xy_2_object()(ps, pt);
      CGAL_assertion(res != EQUAL);
      is_degen = false;
      is_right = (res == SMALLER);

      l = kernel.construct_line_2_object()(ray);
      is_vert = kernel.is_vertical_2_object()(ray);
      is_horiz = kernel.is_horizontal_2_object()(ray);
      has_pos_slope = _has_positive_slope();
    }

    /*! Constructor from a line.
     * \param ln The line.
     * \pre The line is not degenerate.
     */
    _Linear_object_cached_2(const Line_2& ln) :
      l(ln),
      has_source(false),
      has_target(false)
    {
      Kernel kernel;

      CGAL_assertion_msg(! kernel.is_degenerate_2_object()(ln),
                         "Cannot construct a degenerate line.");

      auto construct_vertex = kernel.construct_point_on_2_object();
      ps = construct_vertex(ln, 0);         // Some point on the line.
      pt = construct_vertex(ln, 1);         // Some point further on the line.

      Comparison_result res = kernel.compare_xy_2_object()(ps, pt);
      CGAL_assertion(res != EQUAL);
      is_degen = false;
      is_right = (res == SMALLER);

      is_vert = kernel.is_vertical_2_object()(ln);
      is_horiz = kernel.is_horizontal_2_object()(ln);
      has_pos_slope = _has_positive_slope();
    }

    /*! Check whether the x-coordinate of the left point is infinite.
     * \return ARR_LEFT_BOUNDARY if the left point is near the boundary;
     *         ARR_INTERIOR if the x-coordinate is finite.
     */
    Arr_parameter_space left_infinite_in_x() const
    {
      if (is_vert || is_degen) return (ARR_INTERIOR);

      return (is_right) ?
        (has_source ? ARR_INTERIOR : ARR_LEFT_BOUNDARY) :
        (has_target ? ARR_INTERIOR : ARR_LEFT_BOUNDARY);
    }

    /*! Check whether the y-coordinate of the left point is infinite.
     * \return ARR_BOTTOM_BOUNDARY if the left point is at y = -oo;
     *         ARR_INTERIOR if the y-coordinate is finite.
     *         ARR_TOP_BOUNDARY if the left point is at y = +oo;
     */
    Arr_parameter_space left_infinite_in_y() const
    {
      if (is_horiz || is_degen) return ARR_INTERIOR;

      if (is_vert) {
        return (is_right) ?
          (has_source ? ARR_INTERIOR : ARR_BOTTOM_BOUNDARY) :
          (has_target ? ARR_INTERIOR : ARR_BOTTOM_BOUNDARY);
      }

      if ((is_right && has_source) || (! is_right && has_target))
        return ARR_INTERIOR;

      return (has_pos_slope ? ARR_BOTTOM_BOUNDARY : ARR_TOP_BOUNDARY);
    }

    /*! Check whether the left point is finite.
     */
    bool has_left() const { return (is_right ? has_source : has_target); }

    /*! Obtain the (lexicographically) left endpoint.
     * \pre The left point is finite.
     */
    const Point_2& left() const
    {
      CGAL_precondition(has_left());
      return (is_right ? ps : pt);
    }

    /*! Set the (lexicographically) left endpoint.
     * \param p The point to set.
     * \pre p lies on the supporting line to the left of the right endpoint.
     */
    void set_left(const Point_2& p,
                  bool CGAL_assertion_code(check_validity) = true)
    {
      CGAL_precondition(! is_degen);

      CGAL_precondition_code(Kernel kernel);
      CGAL_precondition
        (Segment_assertions::_assert_is_point_on(p, l, Has_exact_division()) &&
         (! check_validity || ! has_right() ||
          kernel.compare_xy_2_object()(p, right()) == SMALLER));

      if (is_right) {
        ps = p;
        has_source = true;
      }
      else {
        pt = p;
        has_target = true;
      }
    }

    /*! Set the (lexicographically) left endpoint as infinite.
     */
    void set_left()
    {
      CGAL_precondition(! is_degen);

      if (is_right) has_source = false;
      else has_target = false;
    }

    /*! Check whether the x-coordinate of the right point is infinite.
     * \return ARR_RIGHT_BOUNDARY if the right point is near the boundary;
     *         ARR_INTERIOR if the x-coordinate is finite.
     */
    Arr_parameter_space right_infinite_in_x() const
    {
      if (is_vert || is_degen) return ARR_INTERIOR;

      return (is_right) ?
        (has_target ? ARR_INTERIOR : ARR_RIGHT_BOUNDARY) :
        (has_source ? ARR_INTERIOR : ARR_RIGHT_BOUNDARY);
    }

    /*! Check whether the y-coordinate of the right point is infinite.
     * \return ARR_BOTTOM_BOUNDARY if the right point is at y = -oo;
     *         ARR_INTERIOR if the y-coordinate is finite.
     *         ARR_TOP_BOUNDARY if the right point is at y = +oo;
     */
    Arr_parameter_space right_infinite_in_y() const
    {
      if (is_horiz || is_degen) return ARR_INTERIOR;

      if (is_vert) {
        return (is_right) ?
          (has_target ? ARR_INTERIOR : ARR_TOP_BOUNDARY) :
          (has_source ? ARR_INTERIOR : ARR_TOP_BOUNDARY);
      }

      if ((is_right && has_target) || (! is_right && has_source))
          return ARR_INTERIOR;

      return (has_pos_slope ? ARR_TOP_BOUNDARY : ARR_BOTTOM_BOUNDARY);
    }

    /*! Check whether the right point is finite.
     */
    bool has_right() const { return (is_right ? has_target : has_source); }

    /*! Obtain the (lexicographically) right endpoint.
     * \pre The right endpoint is finite.
     */
    const Point_2& right() const
    {
      CGAL_precondition(has_right());
      return (is_right ? pt : ps);
    }

    /*! Set the (lexicographically) right endpoint.
     * \param p The point to set.
     * \pre p lies on the supporting line to the right of the left endpoint.
     */
    void set_right(const Point_2& p,
                   bool CGAL_assertion_code(check_validity) = true)
    {
      CGAL_precondition(! is_degen);
      CGAL_precondition_code(Kernel kernel);
      CGAL_precondition
        (Segment_assertions::_assert_is_point_on(p, l, Has_exact_division()) &&
         (! check_validity || ! has_left() ||
          kernel.compare_xy_2_object()(p, left()) == LARGER));

      if (is_right) {
        pt = p;
        has_target = true;
      }
      else {
        ps = p;
        has_source = true;
      }
    }

    /*! Set the (lexicographically) right endpoint as infinite.
     */
    void set_right()
    {
      CGAL_precondition (! is_degen);

      if (is_right) has_target = false;
      else has_source = false;
    }

    /*! Obtain the supporting line.
     */
    const Line_2& supp_line() const
    {
      CGAL_precondition(! is_degen);
      return (l);
    }

    /*! Check whether the curve is vertical.
     */
    bool is_vertical() const
    {
      CGAL_precondition(! is_degen);
      return (is_vert);
    }

    /*! Check whether the curve is degenerate.
     */
    bool is_degenerate() const { return (is_degen); }

    /*! Check whether the curve is directed lexicographic from left to right
     */
    bool is_directed_right() const { return (is_right); }

    /*! Check whether the given point is in the x-range of the object.
     * \param p The query point.
     * \return (true) is in the x-range of the segment; (false) if it is not.
     */
    bool is_in_x_range(const Point_2& p) const
    {
      Kernel kernel;
      typename Kernel_::Compare_x_2 compare_x = kernel.compare_x_2_object();
      Comparison_result res1;

      if (left_infinite_in_x() == ARR_INTERIOR) {
        // Compare with some point on the curve.
        if (left_infinite_in_y() != ARR_INTERIOR) res1 = compare_x(p, ps);
        else res1 = compare_x(p, left());
      }
      else {
        // p is obviously to the right.
        res1 = LARGER;
      }

      if (res1 == SMALLER) return false;
      else if (res1 == EQUAL) return true;

      Comparison_result res2;

      if (right_infinite_in_x() == ARR_INTERIOR) {
        // Compare with some point on the curve.
        if (right_infinite_in_y() != ARR_INTERIOR) res2 = compare_x(p, ps);
        else res2 = compare_x(p, right());
      }
      else {
        // p is obviously to the right.
        res2 = SMALLER;
      }

      return (res2 != LARGER);
    }

    /*! Check whether the given point is in the y-range of the object.
     * \param p The query point.
     * \pre The object is vertical.
     * \return (true) is in the y-range of the segment; (false) if it is not.
     */
    bool is_in_y_range(const Point_2& p) const
    {
      CGAL_precondition(is_vertical());

      Kernel kernel;
      typename Kernel_::Compare_y_2 compare_y = kernel.compare_y_2_object();
      Arr_parameter_space inf = left_infinite_in_y();
      Comparison_result res1;

      CGAL_assertion(inf != ARR_TOP_BOUNDARY);
      if (inf == ARR_INTERIOR) res1 = compare_y (p, left());
      else res1 = LARGER;           // p is obviously above.

      if (res1 == SMALLER) return false;
      else if (res1 == EQUAL) return true;

      Comparison_result res2;

      inf = right_infinite_in_y();
      CGAL_assertion(inf != ARR_BOTTOM_BOUNDARY);
      if (inf == ARR_INTERIOR) res2 = compare_y(p, right());
      else res2 = SMALLER;          // p is obviously below.

      return (res2 != LARGER);
    }

  private:
    /*! Determine if the supporting line has a positive slope.
     */
    bool _has_positive_slope() const
    {
      if (is_vert) return true;
      if (is_horiz) return false;

      // Construct a horizontal line and compare its slope the that of l.
      Kernel kernel;
      Line_2 l_horiz =
        kernel.construct_line_2_object()(Point_2(0, 0), Point_2(1, 0));
      return (kernel.compare_slope_2_object()(l, l_horiz) == LARGER);
    }
  };

public:
  // Traits objects
  typedef typename Kernel::Point_2              Point_2;
  typedef Arr_linear_object_2<Kernel>           X_monotone_curve_2;
  typedef Arr_linear_object_2<Kernel>           Curve_2;
  typedef unsigned int                          Multiplicity;

public:
  /*! Default constructor.
   */
  Arr_linear_traits_2() {}

  /// \name Basic functor definitions.
  //@{

  /*! A functor that compares the x-coordinates of two points */
  class Compare_x_2 {
  protected:
    typedef Arr_linear_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_2(const Traits& traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_linear_traits_2<Kernel>;

  public:
    /*! Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) < x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      const Kernel& kernel = m_traits;
      return (kernel.compare_x_2_object()(p1, p2));
    }
  };

  /*! Obtain a Compare_x_2 functor. */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(*this); }

  /*! A functor that compares the he endpoints of an $x$-monotone curve. */
  class Compare_endpoints_xy_2{
  public:
    /*! Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param cv The curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
     Comparison_result operator()(const X_monotone_curve_2& xcv) const
    { return (xcv.is_directed_right()) ? (SMALLER) : (LARGER); }
  };

  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  class Trim_2 {
  protected:
    typedef Arr_linear_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Trim_2(const Traits& traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_linear_traits_2<Kernel>;

  public:
    X_monotone_curve_2 operator()(const X_monotone_curve_2 xcv,
                                  const Point_2 src,
                                  const Point_2 tgt)
    {
      /*
       * "Line_segment, line, and ray" will become line segments
       * when trimmed.
       */
      Equal_2 equal = Equal_2();
      Compare_y_at_x_2 compare_y_at_x = m_traits.compare_y_at_x_2_object();

      //preconditions
      //check if source and taget are distinct points and they lie on the line.
      CGAL_precondition(!equal(src, tgt));
      CGAL_precondition(compare_y_at_x(src, xcv) == EQUAL);
      CGAL_precondition(compare_y_at_x(tgt, xcv) == EQUAL);

      //create trimmed line_segment
      X_monotone_curve_2 trimmed_segment;

      if (xcv.is_directed_right() && (tgt.x() < src.x()))
        trimmed_segment = Segment_2(tgt, src);
      else if (! xcv.is_directed_right() && (tgt.x() > src.x()))
        trimmed_segment = Segment_2(tgt, src);
      else trimmed_segment = Segment_2(src, tgt);

      return trimmed_segment;
    }
  };

  Trim_2 trim_2_object() const { return Trim_2(*this); }

  class Construct_opposite_2{
  protected:
    typedef Arr_linear_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Construct_opposite_2(const Traits& traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_linear_traits_2<Kernel>;

  public:
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv) const
    {
      CGAL_precondition(! xcv.is_degenerate());

      X_monotone_curve_2 opp_xcv;

      if (xcv.is_segment()) opp_xcv = Segment_2(xcv.target(), xcv.source());
      if (xcv.is_line()) opp_xcv = Line_2(xcv.get_pt(), xcv.get_ps());
      if (xcv.is_ray()) {
        Point_2 opp_tgt = Point_2( -(xcv.get_pt().x()), -(xcv.get_pt().y()));
        opp_xcv = Ray_2( xcv.source(),  opp_tgt);
      }

      return opp_xcv;
    }
  };

  /*! Get a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(*this); }

  /*! A functor that compares the x-coordinates of two points */
  class Compare_xy_2 {
  public:
    /*! Compare two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2);
     *         SMALLER if x(p1) < x(p2), or if x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL if the two points are equal.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      Kernel kernel;
      return (kernel.compare_xy_2_object()(p1, p2));
    }
  };

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(); }

  /*! A functor that obtains the left endpoint of a segment or a ray. */
  class Construct_min_vertex_2 {
  public:
    /*! Obtain the left endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \pre The left end of cv is a valid (bounded) point.
     * \return The left endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& cv) const
    {
      CGAL_precondition(! cv.is_degenerate());
      CGAL_precondition(cv.has_left());

      return (cv.left());
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  /*! A functor that obtains the right endpoint of a segment or a ray. */
  class Construct_max_vertex_2 {
  public:
    /*! Obtain the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \pre The right end of cv is a valid (bounded) point.
     * \return The right endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& cv) const
    {
      CGAL_precondition(! cv.is_degenerate());
      CGAL_precondition(cv.has_right());

      return (cv.right());
    }
  };

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  /*! A functor that checks whether a given linear curve is vertical. */
  class Is_vertical_2 {
  public:
    /*! Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv) const
    {
      CGAL_precondition(! cv.is_degenerate());
      return (cv.is_vertical());
    }
  };

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(); }

  /*! A functor that compares the y-coordinates of a point and a line at
   * the point x-coordinate
   */
  class Compare_y_at_x_2 {
  protected:
    typedef Arr_linear_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_at_x_2(const Traits& traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_linear_traits_2<Kernel>;

  public:
    /*! Obtain the location of the given point with respect to the input curve.
     * \param cv The curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& cv) const
    {
      CGAL_precondition(! cv.is_degenerate());
      CGAL_precondition(cv.is_in_x_range(p));

      const Kernel& kernel = m_traits;
      if (! cv.is_vertical())
        // Compare p with the segment's supporting line.
        return (kernel.compare_y_at_x_2_object()(p, cv.supp_line()));

      // Compare with the vertical segment's end-points.
      typename Kernel::Compare_y_2 compare_y = kernel.compare_y_2_object();
      const Comparison_result res1 =
        cv.has_left() ? compare_y(p, cv.left()) : LARGER;
      const Comparison_result res2 =
        cv.has_right() ? compare_y(p, cv.right()) : SMALLER;

      return (res1 == res2) ? res1 : EQUAL;
    }
  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(*this); }

  /*! A functor that compares compares the y-coordinates of two linear
   * curves immediately to the left of their intersection point.
   */
  class Compare_y_at_x_left_2 {
  public:
    /*! Compare the y value of two x-monotone curves immediately to the left
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
                                 const Point_2& CGAL_precondition_code(p)) const
    {
      CGAL_precondition(! cv1.is_degenerate());
      CGAL_precondition(! cv2.is_degenerate());

      Kernel                        kernel;

      // Make sure that p lies on both curves, and that both are defined to its
      // left (so their left endpoint is lexicographically smaller than p).
      CGAL_precondition_code(auto compare_xy = kernel.compare_xy_2_object());

      CGAL_precondition
        (Segment_assertions::_assert_is_point_on(p, cv1,
                                                 Has_exact_division()) &&
         Segment_assertions::_assert_is_point_on(p, cv2, Has_exact_division()));

      CGAL_precondition((! cv1.has_left() ||
                         compare_xy(cv1.left(), p) == SMALLER) &&
                        (! cv2.has_left() ||
                         compare_xy(cv2.left(), p) == SMALLER));

      // Compare the slopes of the two segments to determine thir relative
      // position immediately to the left of q.
      // Notice we use the supporting lines in order to compare the slopes,
      // and that we swap the order of the curves in order to obtain the
      // correct result to the left of p.
      return (kernel.compare_slope_2_object()(cv2.supp_line(), cv1.supp_line()));
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(); }

  /*! A functor that compares compares the y-coordinates of two linear
   * curves immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2 {
  public:
    /*! Compare the y value of two x-monotone curves immediately to the right
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
                                 const Point_2& CGAL_precondition_code(p)) const
    {
      CGAL_precondition(! cv1.is_degenerate());
      CGAL_precondition(! cv2.is_degenerate());

      Kernel kernel;

      // Make sure that p lies on both curves, and that both are defined to its
      // right (so their right endpoint is lexicographically larger than p).
      CGAL_precondition_code(auto compare_xy = kernel.compare_xy_2_object());

      CGAL_precondition
        (Segment_assertions::_assert_is_point_on(p, cv1,
                                                 Has_exact_division()) &&
         Segment_assertions::_assert_is_point_on(p, cv2, Has_exact_division()));

      CGAL_precondition((! cv1.has_right() ||
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

  /*! Obtain a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(); }

  /*! A functor that checks whether two points and two linear curves are
   * identical.
   */
  class Equal_2 {
  public:
    /*! Check whether the two x-monotone curves are the same (have the same
     * graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      CGAL_precondition(! cv1.is_degenerate());
      CGAL_precondition(! cv2.is_degenerate());

      Kernel kernel;
      typename Kernel::Equal_2 equal = kernel.equal_2_object();

      // Check that the two supporting lines are the same.
      if (! equal(cv1.supp_line(), cv2.supp_line()) &&
          ! equal(cv1.supp_line(),
                  kernel.construct_opposite_line_2_object()(cv2.supp_line())))
      {
        return false;
      }

      // Check that either the two left endpoints are at infinity, or they
      // are bounded and equal.
      if ((cv1.has_left() != cv2.has_left()) ||
          (cv1.has_left() && ! equal(cv1.left(), cv2.left())))
      {
        return false;
      }

      // Check that either the two right endpoints are at infinity, or they
      // are bounded and equal.
      return ((cv1.has_right() == cv2.has_right()) &&
              (! cv1.has_right() || equal (cv1.right(), cv2.right())));
    }

    /*! Check whether the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    {
      Kernel kernel;
      return (kernel.equal_2_object()(p1, p2));
    }
  };

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object() const { return Equal_2(); }
  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2 {
  public:
    /*! Obtains the parameter space at the end of a line along the x-axis.
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_LEFT_BOUNDARY  - the line approaches the identification arc from
     *                        the right at the line left end.
     *   ARR_INTERIOR       - the line does not approache the identification arc.
     *   ARR_RIGHT_BOUNDARY - the line approaches the identification arc from
     *                        the left at the line right end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                                   Arr_curve_end ce) const
    {
      CGAL_precondition(! xcv.is_degenerate());
      return (ce == ARR_MIN_END) ?
        xcv.left_infinite_in_x() : xcv.right_infinite_in_x();
    }

    /*! Obtains the parameter space at a point along the x-axis.
     * \param p the point.
     * \return the parameter space at p.
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
    /*! Obtains the parameter space at the end of a line along the y-axis .
     * Note that if the line end coincides with a pole, then unless the line
     * coincides with the identification arc, the line end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the line coincides with the identification arc, it is assumed to
     * be smaller than any other object.
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_BOTTOM_BOUNDARY  - the line approaches the south pole at the line
     *                          left end.
     *   ARR_INTERIOR         - the line does not approache a contraction point.
     *   ARR_TOP_BOUNDARY     - the line approaches the north pole at the line
     *                          right end.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2 & xcv,
                                   Arr_curve_end ce) const
    {
      CGAL_precondition (! xcv.is_degenerate());

      return (ce == ARR_MIN_END) ?
        xcv.left_infinite_in_y() : xcv.right_infinite_in_y();
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    Arr_parameter_space operator()(const Point_2 ) const
    { return ARR_INTERIOR; }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }

  /*! A function object that compares the x-limits of arc ends on the
   * boundary of the parameter space
   */
  class Compare_x_at_limit_2 {
  protected:
    typedef Arr_linear_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_x_at_limit_2(const Traits& traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_linear_traits_2<Kernel>;

  public:
    /*! Compare the x-limit of a vertical line at a point with the x-limit of
     * a line end on the boundary at y = +/- oo.
     * \param p the point direction.
     * \param xcv the line, the endpoint of which is compared.
     * \param ce the line-end indicator -
     *            ARR_MIN_END - the minimal end of xc or
     *            ARR_MAX_END - the maximal end of xc.
     * \return the comparison result:
     *         SMALLER - x(p) < x(xc, ce);
     *         EQUAL   - x(p) = x(xc, ce);
     *         LARGER  - x(p) > x(xc, ce).
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the line xcv lies on a boundary, implying
     *      that xcv1 is vertical.
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xcv,
                                 Arr_curve_end ) const
    {
      CGAL_precondition(! xcv.is_degenerate());
      CGAL_precondition(xcv.is_vertical());

      const Kernel& kernel = m_traits;
      return (kernel.compare_x_at_y_2_object()(p, xcv.supp_line()));
    }

    /*! Compare the x-limits of 2 arcs ends on the boundary of the
     * parameter space at y = +/- oo.
     * \param xcv1 the first arc.
     * \param ce1 the first arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second arc.
     * \param ce2 the second arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce1 end of the line xcv1 lies on a boundary, implying
     *      that xcv1 is vertical.
     * \pre the ce2 end of the line xcv2 lies on a boundary, implying
     *      that xcv2 is vertical.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 Arr_curve_end /* ce1 */,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end /*! ce2 */) const
    {
      CGAL_precondition(! xcv1.is_degenerate());
      CGAL_precondition(! xcv2.is_degenerate());
      CGAL_precondition(xcv1.is_vertical());
      CGAL_precondition(xcv2.is_vertical());

      const Kernel& kernel = m_traits;
      const Point_2 p = kernel.construct_point_2_object()(ORIGIN);
      return (kernel.compare_x_at_y_2_object()(p, xcv1.supp_line(),
                                               xcv2.supp_line()));
    }
  };

  /*! Obtain a Compare_x_at_limit_2 function object */
  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(*this); }

  /*! A function object that compares the x-coordinates of arc ends near the
   * boundary of the parameter space
   */
  class Compare_x_near_limit_2 {
  public:
    /*! Compare the x-coordinates of 2 arcs ends near the boundary of the
     * parameter space at y = +/- oo.
     * \param xcv1 the first arc.
     * \param ce1 the first arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second arc.
     * \param ce2 the second arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce end of the line xcv1 lies on a boundary, implying
     *      that xcv1 is vertical.
     * \pre the ce end of the line xcv2 lies on a boundary, implying
     *      that xcv2 is vertical.
     * \pre the the $x$-coordinates of xcv1 and xcv2 at their ce ends are
     *      equal, implying that the curves overlap!
     */
    Comparison_result
    operator()(const X_monotone_curve_2& CGAL_precondition_code(xcv1),
               const X_monotone_curve_2& CGAL_precondition_code(xcv2),
               Arr_curve_end /*! ce2 */) const
    {
      CGAL_precondition(! xcv1.is_degenerate());
      CGAL_precondition(! xcv2.is_degenerate());
      CGAL_precondition(xcv1.is_vertical());
      CGAL_precondition(xcv2.is_vertical());
      return EQUAL;
    }
  };

  /*! Obtain a Compare_x_near_limit_2 function object */
  Compare_x_near_limit_2 compare_x_near_limit_2_object() const
  { return Compare_x_near_limit_2(); }


  /*! A function object that compares the y-limits of arc ends on the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  protected:
    typedef Arr_linear_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     * The constructor is declared private to allow only the functor
     * obtaining function, which is a member of the nesting class,
     * constructing it.
     */
    Compare_y_near_boundary_2(const Traits& traits) : m_traits(traits) {}

    //! Allow its functor obtaining function calling the private constructor.
    friend class Arr_linear_traits_2<Kernel>;

  public:
    /*! Compare the y-limits of 2 lines at their ends on the boundary
     * of the parameter space at x = +/- oo.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the line end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the lines xcv1 and xcv2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      // Make sure both curves are defined at x = -oo (or at x = +oo).
      CGAL_precondition(! xcv1.is_degenerate());
      CGAL_precondition(! xcv2.is_degenerate());
      CGAL_precondition((ce == ARR_MIN_END &&
                         xcv1.left_infinite_in_x() == ARR_LEFT_BOUNDARY &&
                         xcv2.left_infinite_in_x() == ARR_LEFT_BOUNDARY) ||
                        (ce == ARR_MAX_END &&
                         xcv1.right_infinite_in_x() == ARR_RIGHT_BOUNDARY &&
                         xcv2.right_infinite_in_x() == ARR_RIGHT_BOUNDARY));

      // Compare the slopes of the two supporting lines.
      const Kernel& kernel = m_traits;
      const Comparison_result res_slopes =
        kernel.compare_slope_2_object()(xcv1.supp_line(), xcv2.supp_line());

      if (res_slopes == EQUAL) {
        // In case the two supporting line are parallel, compare their
        // relative position at x = 0, which is the same as their position
        // at infinity.
        const Point_2 p = kernel.construct_point_2_object()(ORIGIN);
        return (kernel.compare_y_at_x_2_object()(p, xcv1.supp_line(),
                                                 xcv2.supp_line()));
      }

      // Flip the slope result if we compare at x = -oo:
      return (ce == ARR_MIN_END) ? CGAL::opposite(res_slopes) : res_slopes;
    }
  };

  /*! Obtain a Compare_y_limit_on_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(*this); }

  //@}

  /// \name Functor definitions for supporting intersections.
  //@{

  class Make_x_monotone_2 {
  public:
    /*! Cut the given curve into x-monotone subcurves and insert them into the
     * given output iterator. As segments are always x_monotone, only one
     * object will be contained in the iterator.
     * \param cv The curve.
     * \param oi an output iterator for the result. Its dereference type is a
     *           variant that wraps a \c Point_2 or an \c X_monotone_curve_2
     *           objects.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    {
      // Wrap the segment with a variant.
      typedef boost::variant<Point_2, X_monotone_curve_2>
        Make_x_monotone_result;
      *oi++ = Make_x_monotone_result(cv);
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
    void operator()(const X_monotone_curve_2& cv, const Point_2& p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      CGAL_precondition (! cv.is_degenerate());

      // Make sure that p lies on the interior of the curve.
      CGAL_precondition_code (
        Kernel kernel;
        typename Kernel::Compare_xy_2 compare_xy = kernel.compare_xy_2_object();
      );

      CGAL_precondition
        (Segment_assertions::_assert_is_point_on(p, cv, Has_exact_division()) &&
         (! cv.has_left() || compare_xy(cv.left(), p) == SMALLER) &&
         (! cv.has_right() || compare_xy(cv.right(), p) == LARGER));

      // Perform the split.
      c1 = cv;
      c1.set_right(p);

      c2 = cv;
      c2.set_left(p);
    }
  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object() const { return Split_2(); }

  class Intersect_2 {
  protected:
    typedef Arr_linear_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Intersect_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_linear_traits_2<Kernel>;

  public:
    /*! Find the intersections of the two given curves and insert them into the
     * given output iterator. As two segments may itersect only once, only a
     * single intersection will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    {
      typedef std::pair<Point_2, Multiplicity>          Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;

      CGAL_precondition(! cv1.is_degenerate());
      CGAL_precondition(! cv2.is_degenerate());

      // Intersect the two supporting lines.
      const Kernel& kernel = m_traits;
      auto res = kernel.intersect_2_object()(cv1.supp_line(), cv2.supp_line());

      // The supporting line are parallel lines and do not intersect:
      if (! res) return oi;

      // Check whether we have a single intersection point.
      const Point_2* ip = boost::get<Point_2>(&*res);
      if (ip != nullptr) {
        // Check whether the intersection point ip lies on both segments.
        const bool ip_on_cv1 = cv1.is_vertical() ?
          cv1.is_in_y_range(*ip) : cv1.is_in_x_range(*ip);

        if (ip_on_cv1) {
          const bool ip_on_cv2 = cv2.is_vertical() ?
            cv2.is_in_y_range(*ip) : cv2.is_in_x_range(*ip);

          if (ip_on_cv2) {
            // Create a pair representing the point with its multiplicity,
            // which is always 1 for line segments.
            Intersection_point ip_mult(*ip, 1);
            *oi++ = Intersection_result(ip_mult);
          }
        }
        return oi;
      }

      // In this case, the two supporting lines overlap.
      // We start with the entire cv1 curve as the overlapping subcurve,
      // then clip it to form the true overlapping curve.
      auto compare_xy = kernel.compare_xy_2_object();
      X_monotone_curve_2 ovlp = cv1;

      if (cv2.has_left()) {
        // If the left endpoint of cv2 is to the right of cv1's left endpoint,
        // clip the overlapping subcurve.
        if (! cv1.has_left()) {
          ovlp.set_left (cv2.left(), false);
        }
        else {
          if (compare_xy(cv1.left(), cv2.left()) == SMALLER)
            ovlp.set_left(cv2.left(), false);
        }
      }

      if (cv2.has_right()) {
        // If the right endpoint of cv2 is to the left of cv1's right endpoint,
        // clip the overlapping subcurve.
        if (! cv1.has_right()) {
          ovlp.set_right(cv2.right(), false);
        }
        else {
          if (compare_xy(cv1.right(), cv2.right()) == LARGER)
            ovlp.set_right(cv2.right(), false);
        }
      }

      // Examine the resulting subcurve.
      Comparison_result cmp_res = SMALLER;

      if (ovlp.has_left() && ovlp.has_right())
        cmp_res = compare_xy(ovlp.left(), ovlp.right());

      if (cmp_res == SMALLER) {
        // We have discovered a true overlapping subcurve:
        *oi++ = Intersection_result(ovlp);
      }
      else if (cmp_res == EQUAL) {
        // The two objects have the same supporting line, but they just share
        // a common endpoint. Thus we have an intersection point, but we leave
        // the multiplicity of this point undefined.
        Intersection_point ip_mult(ovlp.left(), 0);
        *oi++ = Intersection_result(ip_mult);
      }

      return oi;
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const { return Intersect_2(*this); }

  class Are_mergeable_2 {
  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      CGAL_precondition(! cv1.is_degenerate());
      CGAL_precondition(! cv2.is_degenerate());

      Kernel kernel;
      typename Kernel::Equal_2 equal = kernel.equal_2_object();

      // Check whether the two curves have the same supporting line.
      if (! equal(cv1.supp_line(), cv2.supp_line()) &&
          ! equal(cv1.supp_line(),
                  kernel.construct_opposite_line_2_object()(cv2.supp_line())))
        return false;

      // Check whether the left endpoint of one curve is the right endpoint of the
      // other.
      return ((cv1.has_right() && cv2.has_left() &&
               equal(cv1.right(), cv2.left())) ||
              (cv2.has_right() && cv1.has_left() &&
               equal(cv2.right(), cv1.left())));
    }
  };

  /*! Obtain an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object () const { return Are_mergeable_2(); }

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2 {
  protected:
    typedef Arr_linear_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_linear_traits_2<Kernel>;

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
      CGAL_precondition(m_traits.are_mergeable_2_object()(cv2, cv1));

      CGAL_precondition(!cv1.is_degenerate());
      CGAL_precondition(!cv2.is_degenerate());

      Equal_2 equal = m_traits.equal_2_object();

      // Check which curve extends to the right of the other.
      if (cv1.has_right() && cv2.has_left() &&
          equal(cv1.right(), cv2.left()))
      {
        // cv2 extends cv1 to the right.
        c = cv1;

        if (cv2.has_right()) c.set_right(cv2.right());
        else c.set_right();      // Unbounded endpoint.
      }
      else {
        CGAL_precondition(cv2.has_right() && cv1.has_left() &&
                          equal(cv2.right(), cv1.left()));

        // cv1 extends cv2 to the right.
        c = cv2;

        if (cv1.has_right()) c.set_right(cv1.right());
        else c.set_right();      // Unbounded endpoint.
      }
    }
  };

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object() const { return Merge_2(*this); }
  //@}

  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  typedef double                          Approximate_number_type;

  class Approximate_2 {
  public:
    /*! Obtain an approximation of a point coordinate.
     * \param p The exact point.
     * \param i The coordinate index (either 0 or 1).
     * \pre i is either 0 or 1.
     * \return An approximation of p's x-coordinate (if i == 0), or an
     *         approximation of p's y-coordinate (if i == 1).
     */
    Approximate_number_type operator()(const Point_2& p, int i) const
    {
      CGAL_precondition((i == 0) || (i == 1));
      return (i == 0) ? CGAL::to_double(p.x()) : CGAL::to_double(p.y());
    }
  };

  /*! Obtain an Approximate_2 functor object. */
  Approximate_2 approximate_2_object() const { return Approximate_2(); }

  class Construct_x_monotone_curve_2 {
  public:
    /*! Obtain an x-monotone curve connecting the two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q must not be the same.
     * \return A segment connecting p and q.
     */
    X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q) const
    {
      Kernel kernel;
      Segment_2 seg = kernel.construct_segment_2_object()(p, q);

      return (X_monotone_curve_2(seg));
    }
  };

  /*! Obtain a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(); }
  //@}
};

/*!
 * \class A representation of a segment, as used by the Arr_segment_traits_2
 * traits-class.
 */
template <typename Kernel_>
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
  /*! Default constructor.
   */
  Arr_linear_object_2() : Base() {}

  /*! Constructor from two points.
   * \param s The source point.
   * \param t The target point.
   * \pre The two points must not be the same.
   */
  Arr_linear_object_2(const Point_2& s, const Point_2& t) : Base(s, t) {}

  /*! Constructor from a segment.
   * \param seg The segment.
   * \pre The segment is not degenerate.
   */
  Arr_linear_object_2(const Segment_2& seg) : Base(seg) {}

  /*! Constructor from a ray.
   * \param ray The segment.
   * \pre The ray is not degenerate.
   */
  Arr_linear_object_2(const Ray_2& ray) : Base(ray) {}

  /*! Constructor from a line.
   * \param line The line.
   * \pre The line is not degenerate.
   */
  Arr_linear_object_2(const Line_2& line) : Base(line) {}

  /*! Check whether the object is actually a segment.
   */
  bool is_segment() const
  { return (! this->is_degen && this->has_source && this->has_target); }

  /*! Cast to a segment.
   * \pre The linear object is really a segment.
   */
  Segment_2 segment() const
  {
    CGAL_precondition(is_segment());

    Kernel kernel;
    Segment_2 seg = kernel.construct_segment_2_object()(this->ps, this->pt);
    return seg;
  }

  /*! Check whether the object is actually a ray.
   */
  bool is_ray() const
  { return (! this->is_degen && (this->has_source != this->has_target)); }

  /*! Cast to a ray.
   * \pre The linear object is really a ray.
   */
  Ray_2 ray() const
  {
    CGAL_precondition(is_ray());

    Kernel kernel;
    Ray_2 ray = (this->has_source) ?
      kernel.construct_ray_2_object()(this->ps, this->l) :
      kernel.construct_ray_2_object()
        (this->pt, kernel.construct_opposite_line_2_object()(this->l));
    return ray;
  }

  /*! Check whether the object is actually a line.
   */
  bool is_line() const
  { return (! this->is_degen && ! this->has_source && ! this->has_target); }

  /*! Cast to a line.
   * \pre The linear object is really a line.
   */
  Line_2 line() const
  {
    CGAL_precondition(is_line());
    return (this->l);
  }

  /*! Get the supporting line.
   * \pre The object is not a point.
   */
  const Line_2& supporting_line() const
  {
    CGAL_precondition(! this->is_degen);
    return (this->l);
  }

  /*!
   * Get the source point.
   * \pre The object is a point, a segment or a ray.
   */
  const Point_2& source() const
  {
    CGAL_precondition(! is_line());

    if (this->is_degen) return (this->ps);      // For a point.
    if (this->has_source) return (this->ps);    // For a segment or a ray.
    else return (this->pt);                     // For a "flipped" ray.
  }

  /*! Get the target point.
   * \pre The object is a point or a segment.
   */
  const Point_2& target() const
  {
    CGAL_precondition(! is_line() && ! is_ray());
    return (this->pt);
  }

  /*! Create a bounding box for the linear object.
   */
  Bbox_2 bbox() const
  {
    CGAL_precondition(this->is_segment());
    Kernel kernel;
    Segment_2 seg = kernel.construct_segment_2_object()(this->ps, this->pt);
    return (kernel.construct_bbox_2_object()(seg));
  }

  // Introducing casting operators instead from a curve to
  // Kernel::Segment_2, Kernel::Ray_2, and Kernel::Line_2 creates an
  // umbiguity. The compiler will barf on the last one, because there are
  // 2 constructors of Kernel::Line_2: one from Kernel::Segment_2 and one
  // from Kernel::Ray_2. Together with the cast to Kernel::Line_2, the
  // compiler will have 3 equivalent options to choose from.
};

/*!
 * Exporter for the segment class used by the traits-class.
 */
template <typename Kernel, typename OutputStream>
OutputStream& operator<<(OutputStream& os,
                         const Arr_linear_object_2<Kernel>& lobj)
{
  // Print a letter identifying the object type, then the object itself.
  if (lobj.is_segment()) os << " S " << lobj.segment();
  else if (lobj.is_ray()) os << " R " << lobj.ray();
  else os << " L " << lobj.line();
  return os;
}

/*! Importer for the segment class used by the traits-class.
 */
template <typename Kernel, typename InputStream>
InputStream& operator>>(InputStream& is, Arr_linear_object_2<Kernel>& lobj)
{
  // Read the object type.
  char c;

  do {
    is >> c;
  } while ((c != 'S' && c != 's') &&
           (c != 'R' && c != 'r') &&
           (c != 'L' && c != 'l'));

  // Read the object accordingly.
  if (c == 'S' || c == 's') {
    typename Kernel::Segment_2  seg;
    is >> seg;
    lobj = seg;
  }
  else if (c == 'R' || c == 'r') {
    typename Kernel::Ray_2 ray;
    is >> ray;
    lobj = ray;
  }
  else {
    typename Kernel::Line_2 line;
    is >> line;
    lobj = line;
  }

  return is;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
