// Copyright (c) 2013 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_FLAT_TORUS_TRAITS_2_H
#define CGAL_ARR_FLAT_TORUS_TRAITS_2_H

/*! \file
 * A class template that handles curves on the flat torus.
 * Any instance of which is suitable as a geometry traits class for the
 * arrangement on surface package.
 */

#include <fstream>

#include <CGAL/config.h>
#include <CGAL/tags.h>
#include <CGAL/intersections.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/number_utils.h>
#include <CGAL/Fraction_traits.h>

namespace CGAL {

#define CGAL_X_MINUS_1_Y_0      0
#define CGAL_X_MINUS_8_Y_6      1
#define CGAL_X_MINUS_11_Y_7     2

#ifndef CGAL_IDENTIFICATION_XY
#define CGAL_IDENTIFICATION_XY  CGAL_X_MINUS_1_Y_0
#endif

template <typename Kernel> class Arr_point_on_flat_torus_3;
template <typename Kernel> class Arr_x_monotone_curve_on_flat_torus_3;
template <typename Kernel> class Arr_curve_on_flat_torus_3;

/*! A traits class-template for constructing and maintaining curves on the
 * flat torus. It is parameterized by a (linear) geometry kernel, which it
 * also derives from.
 */
template <typename Kernel_>
class Arr_flat_torus_traits_2 : public Kernel_ {
public:
  typedef Kernel_                               Kernel;

  // Category tags:
  typedef Tag_true                              Has_left_category;
  typedef Tag_true                              Has_merge_category;
  typedef Tag_false                             Has_do_intersect_category;

  typedef Arr_identified_side_tag               Left_side_category;
  typedef Arr_identified_side_tag               Bottom_side_category;
  typedef Arr_identified_side_tag               Top_side_category;
  typedef Arr_identified_side_tag               Right_side_category;

protected:
  friend class Arr_point_on_flat_torus_3<Kernel_>;
  friend class Arr_x_monotone_curve_on_flat_torus_3<Kernel_>;
  friend class Arr_curve_on_flat_torus_3<Kernel_>;

public:
  /*! Default constructor */
  Arr_flat_torus_traits_2(){}

protected:
  typedef typename Kernel::FT                   FT;

public:
  // Traits objects
  typedef Arr_point_on_flat_torus_3<Kernel>             Point_2;
  typedef Arr_x_monotone_curve_on_flat_torus_3<Kernel>  X_monotone_curve_2;
  typedef Arr_curve_on_flat_torus_3<Kernel>             Curve_2;
  typedef size_t                                        Multiplicity;

public:
  //! \name Basic functor definitions
  //! @{

  //! A functor that compares the x-coordinates of two points.
  class Compare_x_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Compare_x_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Compare the x-coordinates of two points.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return SMALLER - x(p1) < x(p2);
     *         EQUAL   - x(p1) = x(p2);
     *         LARGER  - x(p1) > x(p2).
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     * Note that in the implementation we only check that the points do not
     * lie on the x-boundary, and we use this functor internally to compare
     * the x-coordinates of two points on the y-boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p1));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p2));

      const FT& x1 = p1.x();
      const FT& x2 = p2.x();
      return (x1 < x2) ? SMALLER : ((x1 == x2) ? EQUAL : LARGER);
    }
  };

  /*! Obtain a Compare_x_2 function object.
   * \return an object of type Compare_x_2.
   */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(*this); }

  //! A functor that compares the y-coordinates of two points.
  class Compare_y_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state)
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Compare_y_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Compare the y-coordinates of two points.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return SMALLER - x(p1) < x(p2);
     *         EQUAL   - x(p1) = x(p2);
     *         LARGER  - x(p1) > x(p2).
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     * Note that in the implementation we only check that the points do not
     * lie on the y-boundary, and we use this functor internally to compare
     * the y-coordinates of two points on the x-boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p1));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p2));

      const FT& y1 = p1.y();
      const FT& y2 = p2.y();
      return (y1 < y2) ? SMALLER : ((y1 == y2) ? EQUAL : LARGER);
    }
  };

  /*! Obtain a Compare_y_2 function object.
   * \return an object of type Compare_y_2.
   */
  Compare_y_2 compare_y_2_object() const { return Compare_y_2(*this); }

  /*! A functor that lexigoraphically compares two points;
   * first by x then by y.
   */
  class Compare_xy_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state)
     */
    Compare_xy_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Compare two points lexigoraphically: by x, then by y.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return SMALLER - x(p1) < x(p2);
     *         SMALLER - x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL   - x(p1) = x(p2) and y(p1) = y(p2);
     *         LARGER  - x(p1) = x(p2) and y(p1) > y(p2);
     *         LARGER  - x(p1) > x(p2).
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p1));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p1));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p2));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p2));

      Comparison_result res = m_traits.compare_x_2_object()(p1, p2);
      if (res == EQUAL) return m_traits.compare_y_2_object()(p1, p2);
      return res;
    }
  };

  /*! Obtain a Compare_xy_2 function object.
   * \return an object of type Compare_xy_2.
   */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(*this); }

  //! A functor that obtain the left endpoint of an x-monotone curve.
  class Construct_min_vertex_2 {
  public:
    /*! Obtain the left endpoint of an x-monotone curve.
     * \param xcv the x-monotone curve.
     * \return the left endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& xcv) const
    { return xcv.left(); }
  };

  /*! Obtain a Construct_min_vertex_2 function object.
   * \return an object of type Construct_min_vertex_2.
   */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  //! A functor that obtain the right endpoint of an x-monotone curve.
  class Construct_max_vertex_2 {
  public:
    /*! Obtain the right endpoint of an x-monotone curve.
     * \param xcv the x-monotone curve.
     * \return the right endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& xcv) const
    { return xcv.right(); }
  };

  /*! Obtain a Construct_max_vertex_2 function object.
   * \return an object of type Construct_max_vertex_2.
   */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  /*! A functor that checks whether an x-monotone curve is vertical. */
  class Is_vertical_2 {
  public:
    /*! Check whether a given x-monotone curve is vertical.
     * \param xcv the x-monotone curve.
     * \return true if the curve is vertical; false otherwise.
     * \pre the curve is not degenerate (consists of a single point).
     */
    bool operator()(const X_monotone_curve_2& xcv) const
    {
      CGAL_precondition(!xcv.is_degenerate());
      return xcv.is_vertical();
    }
  };

  /*! Obtain an Is_vertical_2 function object.
   * \return an object of type Is_vertical_2.
   */
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(); }

  /*! A functor that compares the y-coordinates of a point and an
   * x-monotone curve at the point x-coordinate.
   */
  class Compare_y_at_x_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Obtain the location of the given point with respect to the input
     * x-monotone curve.
     * \param xcv the x-monotone curve.
     * \param p the point.
     * \return SMALLER - y(p) < xcv(x(p)), i.e. the point is below the curve;
     *         EQUAL   - p lies on the curve.
     *         LARGER  - y(p) > xcv(x(p)), i.e. the point is above the curve;
     * \pre p lies in the interior of the parameter space.
     * \pre xcv does not lie on the boundary of the parameter space.
     * \pre p lies in the x-range of xcv.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(xcv));
      CGAL_precondition(m_traits.is_in_x_range_object()(p, xcv));

      if (xcv.is_vertical()) {
        // Compare the point with the bottom endpoint. If smaller or equal,
        // return SMALLER or EQUAL, respectively.
        // Otherwise, compare with the top endpoint. If larger, return LARGER.
        // Otherwise, return EQUAL:
        const Point_2& bottom = xcv.bottom();
        const Point_2& top = xcv.bottom();
        if (! bottom.is_on_y_boundary()) {
          Comparison_result res = m_traits.compare_y_2_object()(p, bottom);
          if (res != LARGER) return res;
        }
        if (top.is_on_y_boundary()) return EQUAL;
        Comparison_result res = m_traits.compare_y_2_object()(p, top);
        return (res == LARGER) ? LARGER : EQUAL;
      }

      // Compute xcv(x(p)):
      const Point_2& source = xcv.source();
      const Point_2& target = xcv.target();

      const FT& p_x = p.x();
      const FT& p_y = p.y();
      const FT& source_x = xcv.source_x();
      const FT& source_y = xcv.source_y();
      const FT& target_x = xcv.target_x();
      const FT& target_y = xcv.target_y();

      typename Kernel::Orientation_2 orient = m_traits.orientation_2_object();
      return (xcv.is_directed_right()) ?
        orient(xcv.source_point(), xcv.target_point(), p) :
        orient(xcv.target_point(), xcv.source_point(), p);
    }
  };

  /*! Obtain a Compare_y_at_x_2 function object.
   * \return an object of type Compare_y_at_x_2.
   */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(*this); }

  /*! A functor that compares the y-coordinates of two x-monotone curves
   * immediately to the left of their intersection point.
   */
  class Compare_y_at_x_left_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Compare_y_at_x_left_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Compare the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param xcv1 the first x-monotone curve.
     * \param xcv2 the second x-monotone curve.
     * \param p the intersection point.
     * \return the relative position of xcv1 with respect to xcv2 immdiately to
     *         the left of p: SMALLER, EQUAL, or LARGER.
     * \pre p lies in the interior of the parameter space.
     * \pre xcv1 and xcv2 lie in the interior of the parameter space.
     * \pre the point p lies on both curves.
     * \pre both curves must be defined to the (lexicographically) left of p.
     * \pre the curves are not degenerate.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 const Point_2& p) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv1));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(xcv1));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv2));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(xcv2));
      CGAL_precondition(!xcv1.is_degenerate());
      CGAL_precondition(!xcv2.is_degenerate());
      CGAL_precondition(p == xcv1.right());
      CGAL_precondition(p == xcv2.right());

      // If Both curves are vertical, they overlap:
      if (xcv1.is_vertical() && xcv2.is_vertical()) return EQUAL;
      if (xcv1.is_vertical()) return SMALLER;
      if (xcv2.is_vertical()) return LARGER;

      // None of the curves is verticel.
      // Compare the y-coord. at the x-coord. of the most right left-endpoint.
      const Point_2& left1 = xcv1.left();
      const Point_2& left2 = xcv2.left();

      // Compare the x-coordinate of the left endpoints of the ciurves.
      Comparison_result res = (left1.is_on_x_boundary()) ?
        ((left2.is_on_x_boundary()) ? EQUAL : SMALLER) :
        ((left2.is_on_x_boundary()) ? LARGER :
         m_traits.compare_x_2_object()(left1, left2));

      Arr_parameter_space left_py2 =
        xcv2.is_directed_right() ?
        xcv2.source_parameter_space_in_y() : xcv2.target_parameter_space_in_y();

      if (res == SMALLER) {
        if (left_py2 == ARR_BOTTOM_BOUNDARY) return LARGER;
        if (left_py2 == ARR_TOP_BOUNDARY) return SMALLER;
        // use left2 and xcv1:
        return opposite(m_traits.compare_y_at_x_2_object()(left2, xcv1));
      }

      Arr_parameter_space left_py1 =
        xcv1.is_directed_right() ?
        xcv1.source_parameter_space_in_y() : xcv1.target_parameter_space_in_y();

      if (res == LARGER) {
        if (left_py1 == ARR_BOTTOM_BOUNDARY) return SMALLER;
        if (left_py1 == ARR_TOP_BOUNDARY) return LARGER;
        // use left1 and xcv2:
        return m_traits.compare_y_at_x_2_object()(left1, xcv2);
      }

      // (res == EQUAL)
      if (left1.is_on_x_boundary()) {
        // Both, left1 and left2, are on the left boundary.
        if ((left_py1 == ARR_INTERIOR) && (left_py2 == ARR_INTERIOR))
          return m_traits.compare_y_2_object()(left1, left2);
        if (left_py1 == left_py2) return EQUAL;
        return (left_py1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
      }

      return m_traits.compare_y_2_object()(left1, left2);
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 function object.
   * \return an object of type Compare_y_at_x_left_2.
   */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(*this); }

  /*! A functor that compares the y-coordinates of two x-monotone curves
   * immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_right_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Compare the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param xcv1 the first x-monotone curve.
     * \param xcv2 the second x-monotone curve.
     * \param p the intersection point.
     * \return the relative position of xcv1 with respect to xcv2 immdiately to
     *         the right of p: SMALLER, EQUAL, or LARGER.
     * \pre p lies in the interior of the parameter space.
     * \pre xcv1 and xcv2 lie in the interior of the parameter space.
     * \pre the point p lies on both curves.
     * \pre both curves must be defined to the (lexicographically) right of p.
     * \pre the curves are not degenerate
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 const Point_2& p) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv1));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(xcv1));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv2));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(xcv2));
      CGAL_precondition(!xcv1.is_degenerate());
      CGAL_precondition(!xcv2.is_degenerate());
      CGAL_precondition(p == xcv1.left());
      CGAL_precondition(p == xcv2.left());

      // If Both curves are vertical, they overlap:
      if (xcv1.is_vertical() && xcv2.is_vertical()) return EQUAL;
      if (xcv1.is_vertical()) return LARGER;
      if (xcv2.is_vertical()) return SMALLER;

      // None of the curves are verticel.
      // Compare the y-coord. at the x-coord of the most left right-endpoint.
      const Point_2& right1 = xcv1.right();
      const Point_2& right2 = xcv2.right();

      Arr_parameter_space right_py2 =
        xcv2.is_directed_right() ?
        xcv2.target_parameter_space_in_y() : xcv2.source_parameter_space_in_y();

      Comparison_result res = (right1.is_on_x_boundary()) ?
        ((right2.is_on_x_boundary()) ? EQUAL : LARGER) :
        ((right2.is_on_x_boundary()) ? SMALLER :
         m_traits.compare_x_2_object()(right1, right2));

      if (res == LARGER) {
        if (right_py2 == ARR_BOTTOM_BOUNDARY) return LARGER;
        if (right_py2 == ARR_TOP_BOUNDARY) return SMALLER;
        // use right2 and xcv1:
        return opposite(m_traits.compare_y_at_x_2_object()(right2, xcv1));
      }

      Arr_parameter_space right_py1 =
        xcv1.is_directed_right() ?
        xcv1.target_parameter_space_in_y() : xcv1.source_parameter_space_in_y();

      if (res == SMALLER) {
        if (right_py1 == ARR_BOTTOM_BOUNDARY) return SMALLER;
        if (right_py1 == ARR_TOP_BOUNDARY) return LARGER;
        // use right1 and xcv2:
        return m_traits.compare_y_at_x_2_object()(right1, xcv2);
      }

      // (res == EQUAL)
      if (right1.is_on_x_boundary()) {
        // Both, left1 and left2, are on the left boundary.
        if ((right_py1 == ARR_INTERIOR) && (right_py2 == ARR_INTERIOR))
          return m_traits.compare_y_2_object()(right1, right2);
        if (right_py1 == right_py2) return EQUAL;
        return (right_py1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
      }

      return m_traits.compare_y_2_object()(right1, right2);
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 function object.
   * \return an object of type Compare_y_at_x_right_2.
   */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(*this); }

  /*! A functor that checks whether two directional points and two x-monotone
   * curves are identical.
   */
  class Equal_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Equal_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Determines whether the two points are the same (have the same image).
     * \param p1 the first point.
     * \param p2 the second point.
     * \return true if the two point are the same; false otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    {
      const FT& x1 = p1.x();
      const FT& y1 = p1.y();
      const FT& x2 = p2.x();
      const FT& y2 = p2.y();

      return (
              ((p1.is_on_x_boundary() && p2.is_on_x_boundary()) && (p1.is_on_y_boundary() && p2.is_on_y_boundary())) ||
              ((p1.is_on_x_boundary() && p2.is_on_x_boundary()) && (y1 == y2)) ||
              (((x1 == x2)) && (p1.is_on_y_boundary() && p2.is_on_y_boundary())) ||
              ((x1 == x2) && (y1 == y2))
              );
    }

    /*! Determines whether the two x-monotone curves are the same (have the
     * same image).
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \return true if the two curves are the same; false otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    {
      /* If both are full in the x and in the y dimensions, they must have
       * either
       * the same vertical directions and the same horizontal directions, or
       * the opposite vertical directions and the opposite horizontal directions
       */
      if (xcv1.is_full() && xcv2.is_full())
        return (((xcv1.is_directed_right() == xcv2.is_directed_right()) &&
                 (xcv1.is_directed_top() == xcv2.is_directed_top())) ||
                ((xcv1.is_directed_right() != xcv2.is_directed_right()) &&
                 (xcv1.is_directed_top() != xcv2.is_directed_top())));

      if (xcv1.is_full_x() != xcv2.is_full_x()) return false;
      if (xcv1.is_full_y() != xcv2.is_full_y()) return false;

      return (operator()(xcv1.left(), xcv2.left()) &&
              operator()(xcv1.right(), xcv2.right()));
    }
  };

  /*! Obtain an Equal_2 function object.
   * \return an object of type Equal_2.
   */
  Equal_2 equal_2_object() const { return Equal_2(*this); }
  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2 {
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state)
     */
    Parameter_space_in_x_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Obtains the parameter space at the end of an x-monotone curve along
     * the x-axis.
     * \param xcv the x-monotone curve.
     * \param ce the curve end indicator:
     *     ARR_MIN_END - the minimal end of xcv or
     *     ARR_MAX_END - the maximal end of xcv
     * \return the parameter space at the ce end of the curve xcv.
     *   ARR_LEFT_BOUNDARY  - the curve reaches the x-boundary from the
     *                        right at the curve left end.
     *   ARR_INTERIOR       - the curve does not reach the x-boundary.
     *   ARR_RIGHT_BOUNDARY - the curve reaches the x boundary from the
     *                        left at the curve right end.
     * \pre xcv does not lie on the x-boundary.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv));
      return (ce == ARR_MIN_END) ?
        (xcv.is_directed_right() ?
         xcv.source_parameter_space_in_x() : xcv.target_parameter_space_in_x()) :
        (xcv.is_directed_right() ?
         xcv.target_parameter_space_in_x() : xcv.source_parameter_space_in_x());
    }

    /*! Obtains the x-parameter space at a point along the x-axis.
     * As a convention, if the point lies on the x-boundary, its x-coordinate
     * is assumed to be smaller than the x-coordinate of a point that does not.
     * \param p the point.
     * \return the parameter space at p.
     * \pre p does not lie on the x-boundary.
     */
    Arr_parameter_space operator()(const Point_2 p) const
    {
      return (m_traits.is_on_x_identification_2_object()(p)) ?
        ARR_LEFT_BOUNDARY : ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_x_2 function object.
   * \return an object of type Parameter_space_in_x_2.
   */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(*this); }

  /*! A function object that obtains the parameter space of a geometric
   * entity along the y-axis
   */
  class Parameter_space_in_y_2 {
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state)
     */
    Parameter_space_in_y_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Obtains the parameter space at the end of an x-monotone curve along
     * the y-axis.
     * If the curve lies on the x boundary, it is assumed to lie on the
     * bottom boundary.
     * \param xcv the x-monotone curve.
     * \param ce the curve end indicator:
     *     ARR_MIN_END - the minimal end of xcv or
     *     ARR_MAX_END - the maximal end of xcv
     * \return the parameter space at the ce end of the curve xcv.
     *   ARR_BOTTOM_BOUNDARY  - the curve reaches the y-boundary from above
     *                          at the curve ce end.
     *   ARR_INTERIOR         - the vurve does not reach the y-boundary.
     *   ARR_TOP_BOUNDARY     - the curve reaches the y-boundary from below
     *                          at the curve ce end.
     * \pre xcv does not lie on the y-boundary.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const
    {
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(xcv));
      return (ce == ARR_MIN_END) ?
        (xcv.is_directed_right() ?
         xcv.source_parameter_space_in_y() : xcv.target_parameter_space_in_y()) :
        (xcv.is_directed_right() ?
         xcv.target_parameter_space_in_y() : xcv.source_parameter_space_in_y());
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * As a convention, if the point lies on the x-boundary, its x-coordinate
     * is assumed to be smaller than the x-coordinate of a point that does not.
     * \param p the point.
     * \return the parameter space at p.
     * \pre p does not lie on the y-boundary.
     */
    Arr_parameter_space operator()(const Point_2 p) const
    {
      return (m_traits.is_on_y_identification_2_object()(p)) ?
        ARR_LEFT_BOUNDARY : ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_y_2 function object.
   * \return an object of type Parameter_space_in_y_2.
   */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(*this); }


  /*! A functor that compares the x-coordinate of curve ends on the boundary of
   * the parameter space with curve ends and points.
   */
  class Compare_x_on_boundary_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Compare_x_on_boundary_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Compare the x-coordinate of a point with the x-coordinate of an
     * x-monotone curve end on the boundary.
     * \param point the point.
     * \param xcv the curve, the endpoint of which is compared.
     * \param ce the curve-end indicator -
     *            ARR_MIN_END - the minimal end of xc or
     *            ARR_MAX_END - the maximal end of xc.
     * \return the comparison result:
     *         SMALLER - x(p) < x(xcv, ce);
     *         EQUAL   - x(p) = x(xcv, ce);
     *         LARGER  - x(p) > x(xcv, ce).
     * \pre p does not lie on the x-boundary of the parameter space.
     * \pre xcv does not lie on the x-boundary of the parameter space.
     * \pre the ce end of the curve xcv lies on the y-boundary of the
     *      parameter space.
     */
    Comparison_result operator()(const Point_2& point,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(point));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv));
      CGAL_precondition(m_traits.parameter_space_in_y_2_object()(xcv, ce) !=
                        ARR_INTERIOR);

      if (((ce == ARR_MIN_END) && !xcv.is_directed_right()) ||
          ((ce == ARR_MAX_END) && xcv.is_directed_right()))
      {
        Arr_parameter_space ps_x = xcv.target_parameter_space_in_x();
        return (ps_x == ARR_LEFT_BOUNDARY) ? LARGER :
          ((ps_x == ARR_RIGHT_BOUNDARY) ? SMALLER :
           (m_traits.compare_x_2_object()(point, xcv.target())));
      }
      Arr_parameter_space ps_x = xcv.source_parameter_space_in_x();
      return (ps_x == ARR_LEFT_BOUNDARY) ? LARGER :
        ((ps_x == ARR_RIGHT_BOUNDARY) ? SMALLER :
         (m_traits.compare_x_2_object()(point, xcv.source())));
    }

    /*! Compare the x-coordinates of 2 curve ends near the boundary of the
     * parameter space.
     * \param xcv1 the first curve.
     * \param ce1 the first curve end indicator -
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second curve.
     * \param ce2 the second curve end indicator -
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce1 end of the curve xcv1 lies on the y-boundary.
     * \pre the ce2 end of the curve xcv2 lies on the y-boundary.
     * \pre xcv1 does not lie in the x-identification curve.
     * \pre xcv2 does not lie in the x- identification curve.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv1));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv2));
      CGAL_precondition(m_traits.parameter_space_in_y_2_object()(xcv1, ce1) !=
                        ARR_INTERIOR);
      CGAL_precondition(m_traits.parameter_space_in_y_2_object()(xcv2, ce2) !=
                        ARR_INTERIOR);

      Arr_parameter_space ps1_x;
      const Point_2* p1 = NULL;
      if (((ce1 == ARR_MIN_END) && !xcv1.is_directed_right()) ||
          ((ce1 == ARR_MAX_END) && xcv1.is_directed_right()))
      {
        ps1_x = xcv1.target_parameter_space_in_x();
        p1 = &(xcv1.target());
      }
      else {
        ps1_x = xcv1.source_parameter_space_in_x();
        p1 = &(xcv1.source());
      }

      Arr_parameter_space ps2_x;
      const Point_2* p2 = NULL;
      if (((ce2 == ARR_MIN_END) && !xcv2.is_directed_right()) ||
            ((ce2 == ARR_MAX_END) && xcv2.is_directed_right()))
      {
        ps2_x = xcv2.target_parameter_space_in_x();
        p2 = &(xcv2.target());
      }
      else {
        ps2_x = xcv2.source_parameter_space_in_x();
        p2 = &(xcv2.source());
      }

      if ((ps1_x == ARR_INTERIOR) && (ps2_x == ARR_INTERIOR))
        return m_traits.compare_x_2_object()(*p1, *p2);
      if (ps1_x == ps2_x) return EQUAL;
      return (ps1_x == ARR_LEFT_BOUNDARY) ? SMALLER :
        ((ps1_x == ARR_RIGHT_BOUNDARY) ? LARGER :
         ((ps2_x == ARR_LEFT_BOUNDARY) ? LARGER : SMALLER));
    }

    /*! Compare the x-coordinate of two given points that lie on the
     * y-boundary.
     * \param p1 the first point.
     * \param p2 the second point.
     * \pre p1 lies in the y-boundary.
     * \pre p2 lies in the y-boundary.
     * \pre p1 does not lie in the x-boundary.
     * \pre p2 does not lie in the x-boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      CGAL_precondition(m_traits.is_on_y_identification_2_object()(p1));
      CGAL_precondition(m_traits.is_on_y_identification_2_object()(p2));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p1));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p2));
      return m_traits.compare_x_2_object()(p1, p2);
    }
  };

  /*! Obtain a Compare_x_on_boundary_2 function object.
   * \return an object of type Compare_x_on_boundary_2.
   */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
  { return Compare_x_on_boundary_2(*this); }


  /*! A functor that compares the x-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_boundary_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Compare_x_near_boundary_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Compare the x-coordinates of 2 x-monotone curve ends near the boundary
     * of the parameter space.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param ce the curve end indicator -
     *            ARR_MIN_END - the minimal end of curves or
     *            ARR_MAX_END - the maximal end of curves.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce) < x(xcv2, ce);
     *         EQUAL   - x(xcv1, ce) = x(xcv2, ce);
     *         LARGER  - x(xcv1, ce) > x(xcv2, ce).
     * \pre the ce end of the curve xcv1 lies on the y-boundary.
     * \pre the ce end of the curve xcv2 lies on the y-boundary.
     * \pre the x-coordinates of xcv1 and xcv2 at their ce end are equal.
     * \pre xcv1 does not lie on the x-boundary.
     * \pre xcv2 does not lie on the x-boundary.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    {
      Parameter_space_in_y_2 py = m_traits.parameter_space_in_y_2_object();
      Arr_parameter_space py1 = py(xcv1, ce);
      Arr_parameter_space py2 = py(xcv2, ce);
      CGAL_precondition(py1 != ARR_INTERIOR);
      CGAL_precondition(py2 != ARR_INTERIOR);

      typename Kernel::Point_2 l1, r1;
      if (xcv1.is_directed_right()) {
        l1 = xcv1.source_point();
        r1 = xcv1.target_point();
      }
      else {
        l1 = xcv1.target_point();
        r1 = xcv1.source_point();
      }
      typename Kernel::Point_2 l2, r2;
      if (xcv2.is_directed_right()) {
        l2 = xcv2.source_point();
        r2 = xcv2.target_point();
      }
      else {
        l2 = xcv2.target_point();
        r2 = xcv2.source_point();
      }

      typename Kernel::Orientation_2 orient = m_traits.orientation_2_object();
      if (py1 == ARR_BOTTOM_BOUNDARY) {
        if (ce == ARR_MAX_END) return orient(l2, r1, l1);
        return orient(r2, l1, r1);
      }
      if (ce == ARR_MAX_END) return orient(l1, r1, l2);
      return orient(r1, l1, r2);
    }
  };

  /*! Obtain a Compare_x_near_boundary_2 function object.
   * \return an object of type Compare_x_near_boundary_2.
   */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  { return Compare_x_near_boundary_2(*this); }


  /*! A functor that compares the y-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Compare_y_near_boundary_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Compare the y-coordinates of 2 x-monotone curves at their ends near
     * the boundary of the parameter space.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param ce the curve end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the curves xcv1 and xcv2 lie either on the left
     *      boundary or on the right boundary of the parameter space.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    {
      CGAL_precondition(!xcv1.is_degenerate());
      CGAL_precondition(!xcv2.is_degenerate());
      CGAL_precondition(m_traits.parameter_space_in_x_2_object()(xcv1, ce) !=
                        ARR_INTERIOR);
      CGAL_precondition(m_traits.parameter_space_in_x_2_object()(xcv2, ce) !=
                        ARR_INTERIOR);

      typename Kernel::Point_2 l1, r1;
      if (xcv1.is_directed_right()) {
        l1 = xcv1.source_point();
        r1 = xcv1.target_point();
      }
      else {
        l1 = xcv1.target_point();
        r1 = xcv1.source_point();
      }
      typename Kernel::Point_2 l2, r2;
      if (xcv2.is_directed_right()) {
        l2 = xcv2.source_point();
        r2 = xcv2.target_point();
      }
      else {
        l2 = xcv2.target_point();
        r2 = xcv2.source_point();
      }

      const Kernel& kernel = m_traits;
      typename Kernel::Orientation_2 orient = kernel.orientation_2_object();
      typename Kernel::Compare_y_2 cmp_y = kernel.compare_y_2_object();
      Arr_parameter_space px1 =
        m_traits.parameter_space_in_x_2_object()(xcv1, ce);
      if (px1 == ARR_LEFT_BOUNDARY) {
        CGAL_assertion(ce == ARR_MIN_END);
        Comparison_result rc = cmp_y(l1, l2);
        if (rc != EQUAL) return rc;
        return orient(r1, l1, r2);
      }
      CGAL_assertion(px1 == ARR_RIGHT_BOUNDARY);
      CGAL_assertion(ce == ARR_MAX_END);
      Comparison_result rc = cmp_y(r1, r2);
      if (rc != EQUAL) return rc;
      return orient(l2, r2, l1);
    }
  };

  /*! Obtain a Compare_y_near_boundary_2 function object.
   * \return an object of type Compare_y_near_boundary_2.
   */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(*this); }

  /*! A functor that compares the y-coordinate of two given points
   * that lie on the y-boundary.
   */
  class Compare_y_on_boundary_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_on_boundary_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Compare the y-coordinate of two given points that lie on the
     * x-boundary of the parameter space.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return SMALLER - p1 is lexicographically smaller than p2;
     *         EQUAL   - p1 and p2 coincides;
     *         LARGER  - p1 is lexicographically larger than p2;
     * \pre p1 lies on the x-boundary.
     * \pre p2 lies on the x-boundary.
     * \pre p1 does not lie on the y-boundary.
     * \pre p2 does not lie on the y-boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      CGAL_precondition(m_traits.is_on_x_identification_2_object()(p1));
      CGAL_precondition(m_traits.is_on_x_identification_2_object()(p2));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p1));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p2));
      return m_traits.compare_y_2_object()(p1, p2);
    }
  };

  /*! Obtain a Compare_y_on_boundary_2 function object.
   * \return an object of type Compare_y_on_boundary_2.
   */
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const
  { return Compare_y_on_boundary_2(*this); }

  /*! A functor that indicates whether a geometric object lies on the
   * y-boundary.
   */
  class Is_on_x_identification_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

  public:
    /*! Determine whether a point lies on the y-boundary.
     * \return a Boolean indicating whether p lies on the y-boundary.
     */
    bool operator()(const Point_2& p) const { return p.is_on_x_boundary(); }

    /*! Determine whether an x-monotone curve lies on the y-boundary.
     * \param xcv the curve.
     * \return a Boolean indicating whether xcv lies on the y-boundary.
     */
    bool operator()(const X_monotone_curve_2& xcv) const
    {
      return (!xcv.is_full_x() &&
              operator()(xcv.source()) && operator()(xcv.target()));
    }
  };

  /*! Obtain a Is_on_x_identification_2 function object.
   * \return an object of type Is_on_x_identification_2.
   */
  Is_on_x_identification_2 is_on_x_identification_2_object() const
  { return Is_on_x_identification_2(); }

  /*! A functor that indicates whether a geometric object lies on the
   * x-boundary.
   */
  class Is_on_y_identification_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Is_on_y_identification_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Determine whether a point lies on the x-boundary.
     * \param p the point.
     * \return a Boolean indicating whether p lies on the x-boundary.
     */
    bool operator()(const Point_2& p) const { return p.is_on_y_boundary(); }

    /*! Determine whether an x-monotone curve lies on the x-boundary.
     * \param xcv the curve.
     * \return a Boolean indicating whether xcv lies on the x-boundary.
     */
    bool operator()(const X_monotone_curve_2& xcv) const
    {
      return (!xcv.is_full_y() &&
              operator()(xcv.source()) && operator()(xcv.target()));
    }
  };

  /*! Obtain a Is_on_y_identification_2 function object.
   * \return an object of type Is_on_y_identification_2.
   */
  Is_on_y_identification_2 is_on_y_identification_2_object() const
  { return Is_on_y_identification_2(*this); }
  //@}

  /// \name Functor definitions for supporting intersections.
  //@{

  /*! A functor that divides an curve into x-monotone curves. That are,
   * curves that do not cross the identification curve.
   */
  class Make_x_monotone_2 {
  public:
    typedef CGAL::Fraction_traits<FT>                   My_traits;
    typedef typename My_traits::Numerator_type          Numerator_type;
    typedef typename My_traits::Denominator_type        Denominator_type;
    typedef typename My_traits::Compose                 Compose;
    typedef typename My_traits::Decompose               Decompose;
    typedef typename CGAL::Coercion_traits<Numerator_type,
                                           Denominator_type>::Type
                                                        Result_type;

  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state)
     */
    Make_x_monotone_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

    //! The reciprocal of the curve slope (dx/dy).
    mutable FT m_slope;

    /*! Cut the given curve into x-monotone subcurves and insert them into the
     * given output iterator.
     * \param source the curve source point.
     * \param target the curve target point.
     * \param oi the output iterator, whose value-type is Object. The output
     *           object is a wrapper of an X_monotone_curve_2.
     * \return the past-the-end iterator.
     * \pre The curve is not degenerate.
     * \pre The curve does not cross an x-grid line.
     */
    template<typename OutputIterator_>
    OutputIterator_ operator()(const Point_2& source, const Point_2& target,
                               OutputIterator_ oi) const
    {
      FT xs = source.x();
      FT ys = source.y();
      FT xt = target.x();
      FT yt = target.y();

      Numerator_type ys_num;
      Denominator_type ys_den ;
      Decompose()(ys, ys_num, ys_den);
      Result_type ys_quotient, ys_reminder;
      CGAL::div_mod(ys_num, ys_den, ys_quotient, ys_reminder);

      Numerator_type yt_num;
      Denominator_type yt_den;
      Decompose()(yt, yt_num, yt_den);
      Result_type yt_quotient, yt_reminder;
      CGAL::div_mod(yt_num, yt_den, yt_quotient, yt_reminder);

      typename Traits::Construct_x_monotone_curve_2 ctr =
        m_traits.construct_x_monotone_curve_2_object();

      ys -= ys_quotient;
      yt_quotient -= ys_quotient;

      // first
      if (yt_quotient >= 1) {
        FT dy = 1 - ys;
        FT x_next_y = m_slope * dy + xs;
        Point_2 from(xs, ys);
        Point_2 to(x_next_y, 1);
        X_monotone_curve_2 xcv = ctr(from, to);
        *oi++ = make_object(xcv);
        xs = x_next_y;
        // \todo The decrementor operator is not supported for the Epec kernel?
        // --yt_quotient;
        yt_quotient = yt_quotient - 1;
      }
      // rest
      Point_2 curr(xs, 0);
      while (yt_quotient >= 1) {
        FT x_next_y = m_slope + xs;
        Point_2 to(x_next_y, 1);
        X_monotone_curve_2 xcv = ctr(curr, to);
        *oi++ = make_object(xcv);
        xs = x_next_y;
        // --yt_quotient;
        // \todo The decrementor operator is not supported for the Epec kernel?
        yt_quotient = yt_quotient - 1;
        curr = Point_2(xs, 0);
      }
      // last
      if (yt_reminder != 0) {
        Point_2 to(xt, Compose()(yt_reminder, yt_den));
        X_monotone_curve_2 xcv = ctr(curr, to);
        *oi++ = make_object(xcv);
      }
      return oi;
    }

  public:
    /*! Cut the given curve into x-monotone subcurves and insert them into the
     * given output iterator.
     * \param cv the curve.
     * \param oi the output iterator, whose value-type is Object. The output
     *           object is a wrapper of either an X_monotone_curve_2, or, in
     *           case the input curve is degenerate, a Point_2 object.
     * \return the past-the-end iterator.
     */
    template<typename OutputIterator_>
    OutputIterator_ operator()(const Curve_2& cv, OutputIterator_ oi) const
    {
      if (cv.is_degenerate()) {
        // The curve is a degenerate point---wrap it with an object:
        *oi++ = make_object(cv.right());
        return oi;
      }

      if (cv.is_x_monotone()) {
        // The curve is monotone---wrap it with an object:
        // *oi++ = make_object(X_monotone_curve_2(c));
        const X_monotone_curve_2* xcv = &cv;
        *oi++ = make_object(*xcv);
        return oi;
      }

      const Point_2& source = cv.source();
      const Point_2& target = cv.target();
      FT xs = source.x();
      FT ys = source.y();
      FT xt = target.x();
      FT yt = target.y();

      Numerator_type xs_num;
      Denominator_type xs_den;
      Decompose()(xs, xs_num, xs_den);
      Result_type xs_quotient, xs_reminder;
      CGAL::div_mod(xs_num, xs_den, xs_quotient, xs_reminder);

      Numerator_type xt_num;
      Denominator_type xt_den;
      Decompose()(xt, xt_num, xt_den);
      Result_type xt_quotient, xt_reminder;
      CGAL::div_mod(xt_num, xt_den, xt_quotient, xt_reminder);

      FT dx = xt - xs;
      FT dy = yt - ys;
      FT slope = dy / dx;
      m_slope = dx / dy;

      xs -= xs_quotient;
      xt_quotient -= xs_quotient;

      // first
     if (xt_quotient >= 1) {
        FT dx = 1 - xs;
        FT y_next_x = slope * dx + ys;
        Point_2 from(xs, ys);
        Point_2 to(1, y_next_x);
        *oi++ = operator()(from, to, oi);
        ys = y_next_x;
        // \todo The decrementor operator is not supported for the Epec kernel?
        // --xt_quotient;
        xt_quotient = xt_quotient - 1;
      }
      // rest
      Point_2 curr(0, ys);
      while (xt_quotient >= 1) {
        FT y_next_x = slope + ys;
        Point_2 to(1, y_next_x);
        *oi++ = operator()(curr, to, oi);
        ys = y_next_x;
        // \todo The decrementor operator is not supported for the Epec kernel?
        // --xt_quotient;
        xt_quotient = xt_quotient - 1;
        curr = Point_2(0, ys);
      }
      // last
      if (xt_reminder != 0) {
        Point_2 to(Compose()(xt_reminder, xt_den), yt);
        *oi++ = operator()(curr, to, oi);
      }
      return oi;
    }
  };

  /*! Obtain a Make_x_monotone_2 function object.
   * \return an object of type Make_x_monotone_2.
   */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(*this); }

  /*! A functor that splits an x-monotone curve at a directional point. */
  class Split_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Split_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Split a given x-monotone curve at a given point into two sub-curves.
     * \param xcv the curve to split
     * \param p the split point.
     * \param xcv1 (output) the left resulting subcurve. p is its right
     * endpoint.
     * \param xcv2 (output) the right resulting subcurve. p is its left
     * endpoint.
     * \pre p lies on xcv but is not one of its endpoints.
     * \pre xcv is not degenerate
     */
    void operator()(const X_monotone_curve_2& xcv, const Point_2& p,
                    X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const
    {
      CGAL_precondition(!xcv.is_degenerate());
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p));
      CGAL_precondition(m_traits.orientation_2_object()(xcv1.source_point(),
                                                        xcv2.target_point(),
                                                        p) == COLLINEAR);
      CGAL_precondition(m_traits.orientation_2_object()(xcv2.source_point(),
                                                        xcv2.target_point(),
                                                        p) == COLLINEAR);
      xcv1 = X_monotone_curve_2(xcv.source(), p,
                                xcv.is_horizontal(), xcv.is_vertical(),
                                xcv.is_directed_right(), xcv.is_directed_top());
      xcv2 = X_monotone_curve_2(p, xcv.target(),
                                xcv.is_horizontal(), xcv.is_vertical(),
                                xcv.is_directed_right(), xcv.is_directed_top());
    }
  };

  /*! Obtain a Split_2 function object.
   * \return an object of type Split_2.
   */
  Split_2 split_2_object() const { return Split_2(*this); }

  /*! A functor that computes intersections between x-monotone curves. */
  class Intersect_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Intersect_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

    /*! \class A visitor that handles the object returned by the function that
     * computes the intersection. The return type can be either a point or
     * a segment (in case of overlapping segments).
     */
    template <typename OutputIterator_>
    struct Intersection_visitor {
      typedef void                                      result_type;
      typedef OutputIterator_                           Output_iterator;

      //! The output iterator.
      Output_iterator& m_oi;

      //! The traits.
      const Traits& m_traits;

      /*! Constructor
       * \param oi The output iterator.
       * \param traits The flat torus traits.
       */
      Intersection_visitor(Output_iterator& oi, const Traits& traits) :
        m_oi(oi),
        m_traits(traits)
      {}

      /*! Handle the case where the intersection is a point.
       * \param p The point of intersection.
       */
      result_type operator()(const typename Kernel::Point_2& p) const
      {
        std::pair<typename Traits::Point_2, Multiplicity> po(Point_2(p), 1);
        *m_oi++ = make_object(po);
      }

      /*! Handle the case where the intersection is a segment.
       * \param s The overlapping segment.
       */
      result_type operator()(const typename Kernel::Segment_2& s) const
      {
        typename Traits::Construct_x_monotone_curve_2 ctr =
          m_traits.construct_x_monotone_curve_2_object();
        X_monotone_curve_2 xcv = ctr(Point_2(s.source()), Point_2(s.target()));
        *m_oi++ = make_object(xcv);
      }
    };

  public:
    /*! Find the intersections of the two given curves and insert them into the
     * given output iterator. As two curves may itersect only once,
     * only a single intersection will be inserted into the output iterator.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param oi the output iterator.
     * \return the past-the-end output iterator.
     * \pre xcv1 and xcv2 are not degenerate
     */
    template <typename OutputIterator_>
    OutputIterator_ operator()(const X_monotone_curve_2& xcv1,
                               const X_monotone_curve_2& xcv2,
                               OutputIterator_ oi) const
    {
      CGAL_precondition(!xcv1.is_degenerate());
      CGAL_precondition(!xcv2.is_degenerate());

      typedef OutputIterator_                           Output_iterator;
      typedef typename Kernel::Segment_2                Segment_2;
      typedef typename Kernel::Construct_segment_2      Construct_segment_2;
      typedef typename Kernel::Intersect_2              Intersect_2;

      Construct_segment_2 ctr_seg = m_traits.construct_segment_2_object();
      Segment_2 seg1 = ctr_seg(xcv1.source_point(), xcv1.target_point());
      Segment_2 seg2 = ctr_seg(xcv2.source_point(), xcv2.target_point());
      // std::cout << "segment 1: " << seg1 << std::endl;
      // std::cout << "segment 2: " << seg2 << std::endl;

      const Kernel& kernel = m_traits;
      typename CGAL::cpp11::result_of<Intersect_2(Segment_2, Segment_2)>::type
        result = kernel.intersect_2_object()(seg1, seg2);
      Intersection_visitor<Output_iterator> visitor(oi, m_traits);
      if (result) boost::apply_visitor(visitor, *result);
      return oi;
    }
  };

  /*! Obtain an Intersect_2 function object. */
  Intersect_2 intersect_2_object() const { return Intersect_2(*this); }

  /*! A functor that tests whether two x-monotone curves can be merged. */
  class Are_mergeable_2 {
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Are_mergeable_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \return true if the two curves are mergeable; false otherwise.
     * Two curves are mergeable if:
     * 1. they have the same underlying lines, and
     * 2. share a common endpoint. This point cannot be on the boundary.
     */
    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    {
      if (xcv1.is_empty() || xcv2.is_empty()) return true;

      if (xcv1.is_full_x() || xcv1.is_full_y() ||
          xcv2.is_full_x() || xcv2.is_full_y())
        return false;

      // Merge only curves the directions of which are the same.
      if (xcv1.is_directed_top() != xcv2.is_directed_top()) return false;
      if (xcv1.is_directed_right() != xcv2.is_directed_right()) return false;
      if (xcv1.is_vertical() != xcv2.is_vertical()) return false;

      // Check whether the curves share an endpoint:
      typename Traits::Equal_2 eq = m_traits.equal_2_object();
      const Point_2& source1 = xcv1.source();
      const Point_2& target1 = xcv1.target();
      const Point_2& source2 = xcv2.source();
      const Point_2& target2 = xcv2.target();
      if ((!eq(target1, source2) || target1.is_on_boundary()) &&
          (!eq(target2, source1) || target2.is_on_boundary()))
        return false;

      // Check whether the undelining lines are the same:
      typedef typename Kernel::Line_2                   Line_2;
      typedef typename Kernel::Construct_line_2         Construct_line_2;
      const Kernel& kernel = m_traits;
      Construct_line_2 ctr_line = kernel.construct_line_2_object();
      Line_2 line1 = ctr_line(xcv1.source_point(), xcv1.target_point());
      Line_2 line2 = ctr_line(xcv2.source_point(), xcv2.target_point());
      return kernel.equal_2_object()(line1, line2);
      return false;
    }
  };

  /*! Obtain an Are_mergeable_2 function object.
   * \return an object of type Are_mergeable_2.
   */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(*this); }

  //! A functor that merges two x-monotone curves into one.
  class Merge_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Merge_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Merge two given x-monotone curves into a single curve (spherical_curve).
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param xcv Output: the merged curve.
     * \pre the two curves are mergeable.
     */
    void operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2,
                    X_monotone_curve_2& xcv) const
    {
      CGAL_precondition(m_traits.are_mergeable_2_object()(xcv1, xcv2) == true);

      if (xcv1.is_degenerate() || xcv1.is_empty()) {
        xcv = xcv2;
        return;
      }

      if (xcv2.is_degenerate() || xcv2.is_empty()) {
        xcv = xcv1;
        return;
      }

      // Find the common endpoint:

      // Construct the merged x-monotone curve:
      // Check whether the curves share an endpoint:
      typename Traits::Equal_2 eq = m_traits.equal_2_object();
      typename Traits::Construct_x_monotone_curve_2 ctr =
        m_traits.construct_x_monotone_curve_2_object();

      const Point_2& target1 = xcv1.target();
      const Point_2& source2 = xcv2.source();
      if (eq(target1, source2) && !target1.is_on_boundary()) {
        xcv = ctr(xcv1.source_point(), xcv2.target_point());
        return;
      }
      const Point_2& target2 = xcv2.target();
      const Point_2& source1 = xcv1.source();
      CGAL_assertion(eq(target2, source1) && !target2.is_on_boundary());
      xcv = ctr(xcv2.source_point(), xcv1.target_point());
    }
  };

  /*! Obtain a Merge_2 function object.
   * \return an object of type Merge_2.
   */
  Merge_2 merge_2_object() const { return Merge_2(*this); }
  //@}

  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  typedef double                          Approximate_number_type;

  class Approximate_2 {
  public:
    /*! Obtain an approximation of a point coordinate.
     * \param p the exact point.
     * \param i the coordinate index (either 0 or 1).
     * \pre i is either 0 or 1.
     * \return an approximation of p's x-coordinate (if i == 0), or an
     *         approximation of p's y-coordinate (if i == 1).
     */
    Approximate_number_type operator()(const Point_2& p, int i) const
    {
      CGAL_precondition(i == 0 || i == 1);
      return (i == 0) ? CGAL::to_double(p.x()) : CGAL::to_double(p.y());
    }
  };

  /*! Obtain an Approximate_2 function object.
   * \return an object of type Approximate_2.
   */
  Approximate_2 approximate_2_object() const { return Approximate_2(); }

  //! Construct an x-monotone curve.
  class Construct_x_monotone_curve_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Construct_x_monotone_curve_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Obtain an x-monotone curve connecting the two given endpoints.
     * \param p1 the first point.
     * \param p2 the second point.
     * \pre p1 and p2 must not be the same.
     * \return an x-monotone torusial curve connecting p1 and p2.
     */
    X_monotone_curve_2 operator()(const Point_2& source, const Point_2& target,
                                  bool is_horizontal,
                                  bool is_vertical,
                                  bool is_directed_right,
                                  bool is_directed_top,
                                  bool is_degenerate = false,
                                  bool is_empty = false,
                                  bool is_full_x = false,
                                  bool is_full_y = false) const
    {
      return X_monotone_curve_2(source, target, is_horizontal, is_vertical,
                                is_directed_right, is_directed_top,
                                is_full_x, is_full_y, is_degenerate, is_empty);
    }

    /*! Obtain an x-monotone curve connecting the two given endpoints.
     * \param p1 the first point.
     * \param p2 the second point.
     * \pre p1 and p2 must not be the same.
     * \return an x-monotone torusial curve connecting p1 and p2.
     */
    X_monotone_curve_2 operator()(const Point_2& source,
                                  const Point_2& target) const
    {
      bool is_full_x(((source.x() == 0) && (target.x() == 1)) ||
                     ((target.x() == 0) && (source.x() == 1)));
      bool is_full_y(((source.y() == 0) && (target.y() == 1)) ||
                     ((target.y() == 0) && (source.y() == 1)));
      CGAL_precondition(is_full_x || is_full_y ||
                        !m_traits.equal_2_object()(source, target));

      bool is_degenerate(false);
      bool is_empty(false);

      bool is_directed_right((source.x() < target.x()) ||
                             ((source.x() == target.x()) &&
                              (source.y() < target.y())));
      bool is_directed_top((source.y() < target.y()) ||
                             ((source.y() == target.y()) &&
                              (source.x() < target.x())));
      bool is_horizontal(source.y() == target.y());
      bool is_vertical(source.x() == target.x());
      return X_monotone_curve_2(source, target, is_horizontal, is_vertical,
                                is_directed_right, is_directed_top,
                                is_full_x, is_full_y, is_degenerate, is_empty);
    }
  };

  /*! Obtain a Construct_x_monotone_curve_2 function object.
   * \return an object of type Construct_x_monotone_curve_2.
   */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(*this); }

  //! Construct a curve.
  class Construct_curve_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    /*! Constructor from traits.
     * \param traits the traits (in case it has state).
     */
    Construct_curve_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Obtain a curve connecting the two given endpoints.
     * \param p1 the first point.
     * \param p2 the second point.
     * \pre p1 and p2 must not be the same.
     * \return a curve connecting p1 and p2.
     */
    Curve_2 operator()(const Point_2& source, const Point_2& target) const
    {
      bool is_degenerate(false);
      bool is_empty(false);
      bool is_full_x(((source.x() == 0) && (target.x() == 1)) ||
                     ((target.x() == 1) && (source.x() == 0)));
      bool is_full_y(((source.y() == 0) && (target.y() == 1)) ||
                     ((target.y() == 1) && (source.y() == 0)));
      bool is_directed_right((source.x() < target.x()) ||
                             ((source.x() == target.x()) &&
                              (source.y() < target.y())));
      bool is_directed_top((source.y() < target.y()) ||
                             ((source.y() == target.y()) &&
                              (source.x() < target.x())));
      bool is_horizontal(source.y() == target.y());
      bool is_vertical(source.x() == target.x());
      CGAL_precondition((source.x() != target.x()) ||
                        (source.y() != target.y()));
      return Curve_2(source, target, is_horizontal, is_vertical,
                     is_directed_right, is_directed_top,
                     is_full_x, is_full_y, is_degenerate, is_empty);
    }
  };

  /*! Obtain a Construct_curve_2 function object.
   * \return an object of type Construct_curve_2.
   */
  Construct_curve_2 construct_curve_2_object() const
  { return Construct_curve_2(*this); }

  //@}

  /// \name Functor definitions for the Boolean set-operation traits.
  //@{
  class Compare_endpoints_xy_2 {
  public:
    /*! Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param xcv the curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv)
    { return (xcv.is_directed_right()) ? SMALLER : LARGER; }
  };

  /*! Obtain a Compare_endpoints_xy_2 function object.
   * \return an object of type Compare_endpoints_xy_2.
   */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  class Construct_opposite_2 {
  public:
    /*! Construct an opposite x-monotone (with swapped source and target).
     * \param xcv the curve.
     * \return the opposite curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv)
    { return xcv.opposite(); }
  };

  /*! Obtain a Construct_opposite_2 function object.
   * \retrun an object of type Construct_opposite_2.
   */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }
  //@}

  /*! A functor that determines whether a point is in the x-range of an
   * x-monotone curve.
   */
  class Is_in_x_range {
  public:
    /*! Determine whether the given point is in the x-range of the given curve.
     * \param point the query point.
     * \param xcv the x-monotone curve.
     * \return true if point lies in the x-range of xc; false otherwise.
     * \pre p does not lie on the x-boundary.
     */
    bool operator()(const Point_2& p, const X_monotone_curve_2& xcv) const
    {
      CGAL_precondition(!p.is_on_x_boundary());

      const FT& px = p.x();
      const FT& left_x =
        xcv.is_directed_right() ? xcv.source_x() : xcv.target_x();
      const FT& right_x =
        xcv.is_directed_right() ? xcv.target_x() : xcv.source_x();
      return ((left_x <= px) && (px <= right_x));
    }
  };

  /*! Obtain an Is_in_x_range function object.
   * \return an object of type Is_in_x_range.
   */
  Is_in_x_range is_in_x_range_object() const { return Is_in_x_range(); }

  /*! Extractor for the flat-torus point used by the traits-class.
   */
  template <typename InputStream_>
  InputStream_& read(InputStream_& is, Point_2& point)
  // friend InputStream_& operator>>(InputStream_& is, Point_2& point)
  {
    typename Kernel::Point_2 p;
    is >> p;
    point = Point_2(p);
    return is;
  }

  /*! Extractor for the flat-torus x-monotone curve used by the traits-class.
   */
  template <typename InputStream_>
  InputStream_& read(InputStream_& is, X_monotone_curve_2& xcv)
  // friend InputStream_& operator>>(InputStream_& is, X_monotone_curve_2& xcv)
  {
    Point_2 source, target;
    read(is, source);
    read(is, target);
    xcv = construct_x_monotone_curve_2_object()(source, target);
    return is;
  }

  /*! Extractor for the flat-torus curve used by the traits-class.
   */
  template <typename InputStream_>
  InputStream_& read(InputStream_& is, Curve_2& cv)
  // friend InputStream_& operator>>(InputStream_& is, Curve_2& xcv)
  {
    Point_2 source, target;
    read(is, source);
    read(is, target);
    cv = construct_curve_2_object()(source, target);
    return is;
  }
};

//! Represent a point on a torus using 2 coordinates in the range [0,1].
template <typename Kernel_>
class Arr_point_on_flat_torus_3 : public Kernel_::Point_2 {
public:
  typedef Kernel_                                       Kernel;

  typedef typename Kernel::FT                           FT;

private:
  /*! A flag that indicates whether the point lie on the x-boundary. */
  bool m_on_x_boundary;

  /*! A flag that indicates whether the point lie on the y-boundary. */
  bool m_on_y_boundary;

public:
  //! Default constructor
  Arr_point_on_flat_torus_3() {}

  /*! Constructor from 2 coordinates.
   * \param x the x-coordinate.
   * \param y the y-coordinate.
   */
  Arr_point_on_flat_torus_3(const FT& x, const FT& y) :
    Kernel::Point_2(x, y),
    m_on_x_boundary((x == 0) || (x == 1)),
    m_on_y_boundary((y == 0) || (y == 1))
  {}

  /*! Constructor from a (kernel) point.
   * \param x the x-coordinate.
   * \param y the y-coordinate.
   */
  Arr_point_on_flat_torus_3(const typename Kernel::Point_2& p) :
    Kernel::Point_2(p),
    m_on_x_boundary((p.x() == 0) || (p.x() == 1)),
    m_on_y_boundary((p.y() == 0) || (p.y() == 1))
  {}

  /*! Determine whether the point lies on the x-boundary of the parameter space.
   * \return a Boolean that indicates whether the point lies on the x-boundary
   * of the parameter space.
   */
  bool is_on_x_boundary() const { return m_on_x_boundary; }

  /*! Determine whether the point lies on the y-boundary of the parameter space.
   * \return a Boolean that indicates whether the point lies on the y-boundary
   * of the parameter space.
   */
  bool is_on_y_boundary() const { return m_on_y_boundary; }

  /*! Determine whether the point lies on the boundary of the parameter space.
   * \return a Boolean that indicates whether the point lies on the boundary
   * of the parameter space.
   */
  bool is_on_boundary() const
  { return is_on_x_boundary() || is_on_y_boundary(); }
};

//! A Representation of an x-monotone great circular curve embedded on a torus,
// as used by the Arr_flat_torus_traits_2 traits-class template.
// An x-monotone great circular curve cannot cross boundary of the parameter
// space.
template <typename Kernel_>
class Arr_x_monotone_curve_on_flat_torus_3 {
public:
  typedef Kernel_                                   Kernel;

  typedef typename Kernel::FT                       FT;

protected:
  // For some reason compilation under Windows fails without the qualifier
  typedef CGAL::Arr_point_on_flat_torus_3<Kernel>   Arr_point_on_flat_torus_3;

  //! The source point of the curve.
  Arr_point_on_flat_torus_3 m_source;

  //! The target point of the curve.
  Arr_point_on_flat_torus_3 m_target;

  //! The curve is horizontal.
  bool m_is_horizontal;

  //! The curve is vertical.
  bool m_is_vertical;

  //! Target is xy-lexicographically larger than source.
  bool m_is_directed_right;

  //! Target is yx-lexicographically larger than source.
  bool m_is_directed_top;

  //! The curve spans the entire parameter space along the x-direction.
  bool m_is_full_x;

  //! The curve spans the entire parameter space along the y-direction.
  bool m_is_full_y;

  // The curve is degenerate---it consists of a single point.
  bool m_is_degenerate;

  //! The curve is empty.
  bool m_is_empty;

public:
  //! Default constructor - constructs an empty curve.
  Arr_x_monotone_curve_on_flat_torus_3() :
    m_is_horizontal(false),
    m_is_vertical(false),
    m_is_directed_right(false),
    m_is_directed_top(false),
    m_is_full_x(false),
    m_is_full_y(false),
    m_is_degenerate(false),
    m_is_empty(true)
  {}

  /*! Constructor
   * \param src the source point of the curve
   * \param trg the target point of the curve
   * \param is_horizontal indicates whether the curve is horizontal.
   * \param is_vertical indicates whether the curve is vertical.
   * \param is_directed_right indicates whether the curve is directed
   *        from left to right?
   * \param is_directed_top indicates whether the curve is directed
   *        from bottom to top?
   * \param is_full indicates whether the curve is full.
   * \param is_degenerate indicates whether the curve is degenerate
   *         (a single point)?
   */
  Arr_x_monotone_curve_on_flat_torus_3(const Arr_point_on_flat_torus_3& src,
                                       const Arr_point_on_flat_torus_3& trg,
                                       bool is_horizontal,
                                       bool is_vertical,
                                       bool is_directed_right,
                                       bool is_directed_top,
                                       bool is_full_x = false,
                                       bool is_full_y = false,
                                       bool is_degenerate = false,
                                       bool is_empty = false) :
    m_source(src),
    m_target(trg),
    m_is_horizontal(is_horizontal),
    m_is_vertical(is_vertical),
    m_is_directed_right(is_directed_right),
    m_is_directed_top(is_directed_top),
    m_is_full_x(is_full_x),
    m_is_full_y(is_full_y),
    m_is_degenerate(is_degenerate),
    m_is_empty(is_empty)
  {}

  /*! Set the source endpoint.
   * \param p the endpoint to set.
   */
  void set_source(const Arr_point_on_flat_torus_3& p) { m_source = p; }

  /*! Set the target endpoint.
   * \param p the endpoint to set.
   */
  void set_target(const Arr_point_on_flat_torus_3& p) { m_target = p; }

  void set_is_horizontal(bool flag) { m_is_horizontal = flag; }
  void set_is_vertical(bool flag) { m_is_vertical = flag; }
  void set_is_directed_right(bool flag) { m_is_directed_right = flag; }
  void set_is_directed_top(bool flag) { m_is_directed_top = flag; }
  void set_is_full_x(bool flag) { m_is_full_x = flag; }
  void set_is_full_y(bool flag) { m_is_full_y = flag; }
  void set_is_degenerate(bool flag) { m_is_degenerate = flag; }
  void set_is_empty(bool flag) { m_is_empty = flag; }

  /*! Obtain the source */
  const Arr_point_on_flat_torus_3& source() const { return m_source; }

  /*! Obtain the target */
  const Arr_point_on_flat_torus_3& target() const { return m_target; }

  /*! Obtain the xy-lexicographically left endpoint direction */
  const Arr_point_on_flat_torus_3& left() const
  { return (m_is_directed_right ? m_source : m_target); }

  /*! Obtain the xy-lexicographically right endpoint */
  const Arr_point_on_flat_torus_3& right() const
  { return (m_is_directed_right ? m_target : m_source); }

  /*! Obtain the yx-lexicographically bottom endpoint direction */
  const Arr_point_on_flat_torus_3& bottom() const
  { return (m_is_directed_top ? m_source : m_target); }

  /*! Obtain the yx-lexicographically top endpoint */
  const Arr_point_on_flat_torus_3& top() const
  { return (m_is_directed_top ? m_target : m_source); }

  /*! Determine whether the curve is horizontal. */
  bool is_horizontal() const { return m_is_horizontal; }

  /*! Determine whether the curve is vertical. */
  bool is_vertical() const { return m_is_vertical; }

  /*! Determine whether the curve is directed xy-lexicographically from left
   * to right.
   */
  bool is_directed_right() const { return m_is_directed_right; }

  /*! Determine whether the curve is directed yx-lexicographically from bottom
   * to top.
   */
  bool is_directed_top() const { return m_is_directed_top; }

  /*! Determine whether the curve spans the entire x-parameter_space. */
  bool is_full_x() const { return m_is_full_x; }

  /*! Determine whether the curve spans the entire y-parameter_space. */
  bool is_full_y() const { return m_is_full_y; }

  /*! Determine whether the curve spans the entire x- and y-parameter_spaces. */
  bool is_full() const { return (is_full_x() && is_full_y()); }

  /*! Determine whether the curve is degenerate */
  bool is_degenerate() const { return m_is_degenerate; }

  /*! Determine whether the curve is degenerate */
  bool is_empty() const { return m_is_empty; }

  /*! Determine whether the curve lies on the boundary. */
  bool is_on_boundary() const
  {
    return (source_parameter_space_in_x() == target_parameter_space_in_x()) ||
      (source_parameter_space_in_y() == target_parameter_space_in_y());
  }

  /*! Flip the curve (swap it source and target) */
  Arr_x_monotone_curve_on_flat_torus_3 opposite() const
  {
    Arr_x_monotone_curve_on_flat_torus_3 opp;
    opp.m_source = this->m_target;
    opp.m_target = this->m_source;
    opp.m_is_directed_right = !(this->is_directed_right());
    opp.m_is_directed_top = !(this->is_directed_top());
    opp.m_is_horizontal = this->is_horizontal();
    opp.m_is_vertical = this->is_vertical();
    opp.m_is_full_x = this->is_full_x();
    opp.m_is_full_y = this->is_full_y();
    opp.m_is_degenerate = this->is_degenerate();
    opp.m_is_empty = this->is_empty();
    return opp;
  }

  /*! Obtain the parameter space in x of the source point.
   */
  Arr_parameter_space source_parameter_space_in_x() const
  {
    return (!m_source.is_on_x_boundary()) ? ARR_INTERIOR :
      (is_directed_right() ? ARR_LEFT_BOUNDARY : ARR_RIGHT_BOUNDARY) ;
  }

  /*! Obtain the parameter space in y of the source point.
   */
  Arr_parameter_space source_parameter_space_in_y() const
  {
    return (!m_source.is_on_y_boundary()) ? ARR_INTERIOR :
      (is_directed_top() ? ARR_BOTTOM_BOUNDARY : ARR_TOP_BOUNDARY) ;
  }

  /*! Obtain the parameter space in x of the target point.
   */
  Arr_parameter_space target_parameter_space_in_x() const
  {
    return (!m_target.is_on_x_boundary()) ? ARR_INTERIOR :
      (is_directed_right() ? ARR_RIGHT_BOUNDARY : ARR_LEFT_BOUNDARY) ;
  }

  /*! Obtain the parameter space in y of the target point.
   */
  Arr_parameter_space target_parameter_space_in_y() const
  {
    return (!m_target.is_on_y_boundary()) ? ARR_INTERIOR :
      (is_directed_top() ? ARR_TOP_BOUNDARY : ARR_BOTTOM_BOUNDARY) ;
  }

  /*! Obtain the x-coordinate of the source point.
   * \return the x-coordinate of the source point.
   * The following four functions cannot return a reference because for some
   * kernels, e.g., the Epec kernel, the x() and y() member functions of the
   * point do not return a reference. (They return a local variable.)
   * \todo Return a reference when possible. Is it possible?
   */
  const FT source_x() const
  {
    static FT zero(0);
    static FT one(1);
    return (m_source.is_on_x_boundary() ?
            (is_vertical() ? zero : ((is_directed_right()) ? zero : one)) :
            m_source.x());
  }

  /*! Obtain the y-coordinate of the source point.
   * \return the y-coordinate of the source point.
   */
  const FT source_y() const
  {
    static FT zero(0);
    static FT one(1);
    return (m_source.is_on_y_boundary() ?
            (is_horizontal() ? zero : (is_directed_top() ? zero : one)) :
            m_source.y());
  }

  /*! Obtain the x-coordinate of the target point.
   * \return the x-coordinate of the target point.
   */
  const FT target_x() const
  {
    static FT zero(0);
    static FT one(1);
    return (m_target.is_on_x_boundary() ?
            (is_vertical() ? zero : (is_directed_right() ? one : zero)) :
            m_target.x());
  }

  /*! Obtain the y-coordinate of the target point.
   * \return the y-coordinate of the target point.
   */
  const FT target_y() const
  {
    static FT zero(0);
    static FT one(1);
    return (m_target.is_on_y_boundary() ?
            (is_horizontal() ? zero : (is_directed_top() ? one : zero)) :
            m_target.y());
  }

  /*! Obtain the source point.
   * \return the source point image
   */
  typename Kernel::Point_2 source_point() const
  { return typename Kernel::Point_2(source_x(), source_y()); }

  /*! Obtain the target point.
   * \return the target point image
   */
  typename Kernel::Point_2 target_point() const
  { return typename Kernel::Point_2(target_x(), target_y()); }

  /*! Obtain the segment
   * \return the segment;
   */
  typename Kernel::Segment_2 segment() const
  { return typename Kernel::Segment_2(source_point, target_point); }
};

//! A representation of a geodesic curve embedded on the flat torus,
// used by the Arr_flat_torus_traits_2 traits-class template.
// An curve is uniqely represented by two endpoints, the source s and the
// target t. The points of the curve are the locus of points visited when
// moving from the source s toward the target t in a straight line in the
// parameter space.
template <typename Kernel_>
class Arr_curve_on_flat_torus_3 :
  public Arr_x_monotone_curve_on_flat_torus_3<Kernel_> {
public:
  typedef Kernel_                                    Kernel;

protected:
  typedef Arr_x_monotone_curve_on_flat_torus_3<Kernel> Base;

  typedef typename Base::Arr_point_on_flat_torus_3   Arr_point_on_flat_torus_3;

public:
  //! Default constructor - constructs an empty curve.
  Arr_curve_on_flat_torus_3() : Base() {}

  /*! Constructor
   * \param src the source point of the curve
   * \param trg the target point of the curve
   * \param is_horizontal indicates whether the curve is horizontal.
   * \param is_vertical indicates whether the curve is vertical.
   * \param is_directed_right indicates whether the curve is directed
   *        from left to right.
   * \param is_directed_top indicates whether the curve directed
   *        from bottom to top.
   * \param is_full_x indicates whether the curve is x-full.
   * \param is_full_y indicates whether the curve is y-full.
   * \param is_degenerate indicates whether the curve is degenerate
   *        (a single point)?
   * \param is_empty indicates whether the curve is empty.
   */
  Arr_curve_on_flat_torus_3(const Arr_point_on_flat_torus_3& src,
                            const Arr_point_on_flat_torus_3& trg,
                            bool is_horizontak,
                            bool is_vertical,
                            bool is_directed_right,
                            bool is_directed_top,
                            bool is_full_x = false,
                            bool is_full_y = false,
                            bool is_degenerate = false,
                            bool is_empty = false) :
    Base(src, trg, is_horizontak, is_vertical,
         is_directed_right, is_directed_top,
         is_full_x, is_full_y, is_degenerate, is_empty)
  {}

  /*! Indicates whether the curve is x-monotone
   * \return true if the curve is x-monotone; false otherwise
   */
  bool is_x_monotone() const
  {
    typedef typename Kernel::FT FT;
    FT xs = this->source().x();
    FT ys = this->source().y();
    FT xt = this->target().x();
    FT yt = this->target().y();
    // std::cout << "(" << xs << "," << ys << ")(" << xt << "," << yt << ")"
    //           << std::endl;

    typedef CGAL::Fraction_traits<FT>                   My_traits;
    typedef typename My_traits::Numerator_type          Numerator_type;
    typedef typename My_traits::Denominator_type        Denominator_type;
    typedef typename My_traits::Decompose               Decompose;

    Numerator_type xs_num, xt_num;
    Numerator_type ys_num, yt_num;
    Denominator_type xs_den, xt_den;
    Denominator_type ys_den, yt_den;
    Decompose()(xs, xs_num, xs_den);
    Decompose()(ys, ys_num, ys_den);
    Decompose()(xt, xt_num, xt_den);
    Decompose()(yt, yt_num, yt_den);
    // std::cout << "xs_num " << xs_num << std::endl;
    // std::cout << "xs_den " << xs_den << std::endl;
    // std::cout << "ys_num " << ys_num << std::endl;
    // std::cout << "ys_den " << ys_den << std::endl;

    typedef typename CGAL::Coercion_traits<Numerator_type,
                                           Denominator_type>::Type
                                                        Result_type;

    Result_type xs_quotient, xs_reminder;
    CGAL::div_mod(xs_num, xs_den, xs_quotient, xs_reminder);

    Result_type xt_quotient, xt_reminder;
    CGAL::div_mod(xt_num, xt_den, xt_quotient, xt_reminder);

    Result_type ys_quotient, ys_reminder;
    CGAL::div_mod(ys_num, ys_den, ys_quotient, ys_reminder);

    Result_type yt_quotient, yt_reminder;
    CGAL::div_mod(yt_num, yt_den, yt_quotient, yt_reminder);

    //! \todo use directed_right
    if (xs_quotient > xt_quotient) {
      std::swap(xs_quotient, xt_quotient);
      std::swap(xs_reminder, xt_reminder);
    }

    if ((xs_quotient + 1) < xt_quotient) return false;
    if ((xs_quotient != xt_quotient) && (xt_reminder != 0)) return false;

    //! \todo use directed_top
    if (ys_quotient > yt_quotient) {
      std::swap(ys_quotient, yt_quotient);
      std::swap(ys_reminder, yt_reminder);
    }
    if ((ys_quotient + 1) < yt_quotient) return false;
    if ((ys_quotient != yt_quotient) && (yt_reminder != 0)) return false;

    return true;
  }
};

/*! Inserter for the flat-torus point used by the traits-class.
 */
template <typename Kernel_, typename InputStream_>
InputStream_& operator>>(InputStream_& os,
                         const Arr_point_on_flat_torus_3<Kernel_>& point)
{
  typedef Kernel_       Kernel;
  const typename Kernel::Point_2* p = &point;
  os << *p;
  return os;
}

/*! Inserter for the flat-torus x-monotone curve used by the traits-class.
 */
template <typename Kernel_, typename OutputStream_>
OutputStream_&
operator<<(OutputStream_& os,
           const Arr_x_monotone_curve_on_flat_torus_3<Kernel_>& xcv)
{
  os << "(";
  os << xcv.source() << ", " << xcv.target();
#if defined(CGAL_ARR_CURVE_ON_FLAT_TORUS_DETAILS)
#endif
  os << ")";
  return os;
}

/*! Inserter for the flat-torus curve used by the traits-class.
 */
template <typename Kernel_, typename OutputStream_>
OutputStream_&
operator<<(OutputStream_& os, const Arr_curve_on_flat_torus_3<Kernel_>& cv)
{
  os << "(";
  os << cv.source() << ", " << cv.target();
#if defined(CGAL_ARR_CURVE_ON_FLAT_TORUS_DETAILS)
#endif
  os << ")";
  return os;
}

} //namespace CGAL

#endif
