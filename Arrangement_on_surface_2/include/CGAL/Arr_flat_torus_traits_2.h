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
 * A class template that handles geodesics arcs embedded on the flat torus.
 * Any instance of which is suitable as a geometry traits class for the
 * arrangement on surface package.
 */

#include <fstream>

#include <CGAL/config.h>
#include <CGAL/tags.h>
#include <CGAL/intersections.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>

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

/*! A traits class-template for constructing and maintaining arcs of great
 * circles embedded on tori. It is parameterized from a (linear) geometry
 * kernel, which it also derives from
 */
template <typename Kernel_>
class Arr_flat_torus_traits_2 : public Kernel_ {
public:
  typedef Kernel_                              Kernel;

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

    //! The traits (in case it has state)
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
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
     * lie on the x-boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p1));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p2));

      const FT& d1 = p1.x();
      const FT& d2 = p2.x();
      return (d1 < d2) ? SMALLER : ((d1 == d2) ? EQUAL : LARGER);
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

    /*! Constructor
     * \param traits the traits (in case it has state)
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
     * lie on the y-boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p1));
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(p2));

      const FT& d1 = p1.y();
      const FT& d2 = p2.y();
      return (d1 < d2) ? SMALLER : ((d1 == d2) ? EQUAL : LARGER);
    }
  };

  /*! Obtain a Compare_y_2 function object.
   * \return an object of type Compare_y_2.
   */
  Compare_x_2 compare_y_2_object() const { return Compare_y_2(*this); }

  /*! A functor that compares two directional points lexigoraphically:
   * by x, then by y.
   */
  class Compare_xy_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
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

      Comparison_result m_traits.compare_x_object()(p1, p2);
      if (res == EQUAL) return m_traits.compare_y_object()(p1, p2);
      return res;
    }
  };

  /*! Obtain a Compare_xy_2 function object.
   * \return an object of type Compare_xy_2.
   */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(*this); }

  /*! A functor that obtain the left endpoint of an x-monotone curve. */
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

  /*! A functor that obtain the right endpoint of an x-monotone curve. */
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

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Return the location of the given point with respect to the input
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
        Comparison_result res = m_traits->compare_y_object()(p, xcv.bottom());
        if (res != LARGER) return cr;
        res = m_traits->compare_y_object()(p, xcv.top());
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
      FT y = source_y +
        (p_y - source_x) * (left_y - target_y) / (source_x - target_x);
      return (p_y < y) ? SMALLER : ((p_y == y) ? EQUAL : LARGER);
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

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
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

      Arr_parameter_space left_py1 =
        xcv1.is_directed_right() ?
        xcv1.source_parameter_space_in_y() : xcv1.target_parameter_space_in_y();

      Arr_parameter_space left_py2 =
        xcv2.is_directed_right() ?
        xcv1.source_parameter_space_in_y() : xcv1.target_parameter_space_in_y();

      Comparison_result res = (left1.is_on_x_boundary()) ?
        ((left2.is_on_x_boundary()) ? EQUAL : SMALLER) :
        ((left2.is_on_x_boundary()) ? LARGER :
         m_traits->compare_x_object()(left1, left2));

      if (res == SMALLER) {
        if (left_py2 == ARR_BOTTOM_BOUNDARY) return LARGER;
        if (left_py2 == ARR_TOP_BOUNDARY) return SMALLER;
        // use left2 and xcv1:
        return m_traits->compare_y_at_x_object()(left2, xcv1);
      }

      if (res == LERGER) {
        if (left_py1 == ARR_BOTTOM_BOUNDARY) return SMALLER;
        if (left_py1 == ARR_TOP_BOUNDARY) return LARGER;
        // use left1 and xcv2:
        return m_traits->compare_y_at_x_object()(left1, xcv2);
      }

      // (res == EQUAL)
      if (left1.is_on_x_boundary()) {
        // Both, left1 and left2, are on the left boundary.
        if ((left_py1 == ARR_INTERIOR) && (left_py2 == ARR_INTERIOR))
          return m_traits->compare_y_object()(left1, left2);
        if (left_py1 == left_py2) return EQUAL;
        return (left_py1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
      }

      return m_traits->compare_y_object()(left1, left2);
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 function object.
   * \return an object of type Compare_y_at_x_left_2.
   */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(*this); }

  //! A functor that compares the y-coordinates of two x-monotone curves
  // immediately to the right of their intersection point.
  class Compare_y_at_x_right_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
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

      // Non of the arc is verticel.
      // Compare the y-coord. at the x-coord of the most left right-endpoint.
      const Point_2& right1 = xcv1.right();
      const Point_2& right2 = xcv2.right();

      Arr_parameter_space right_py1 =
        xcv1.is_directed_right() ?
        xcv1.target_parameter_space_in_y() : xcv1.source_parameter_space_in_y();

      Arr_parameter_space right_py2 =
        xcv2.is_directed_right() ?
        xcv1.target_parameter_space_in_y() : xcv1.source_parameter_space_in_y();

      Comparison_result res = (right1.is_on_x_boundary()) ?
        ((right2.is_on_x_boundary()) ? EQUAL : LARGER) :
        ((right2.is_on_x_boundary()) ? SMALLER :
         m_traits->compare_x_object()(right1, right2));

      if (res == LARGER) {
        if (right_py2 == ARR_BOTTOM_BOUNDARY) return LARGER;
        if (right_py2 == ARR_TOP_BOUNDARY) return SMALLER;
        // use right2 and xcv1:
        return m_traits->compare_y_at_x_object()(right2, xcv1);
      }

      if (res == SMALLER) {
        if (right_py1 == ARR_BOTTOM_BOUNDARY) return SMALLER;
        if (right_py1 == ARR_TOP_BOUNDARY) return LARGER;
        // use left1 and xcv2:
        return m_traits->compare_y_at_x_object()(left1, xcv2);
      }

      // (res == EQUAL)
      if (right1.is_on_x_boundary()) {
        // Both, left1 and left2, are on the left boundary.
        if ((right_py1 == ARR_INTERIOR) && (right_py2 == ARR_INTERIOR))
          return m_traits->compare_y_object()(left1, left2);
        if (right_py1 == right_py2) return EQUAL;
        return (right_py1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
      }

      return m_traits->compare_y_object()(left1, left2);
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

    /*! The traits (in case it has state) */
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

      return (((p1.is_on_x_boundary() && p2.is_on_x_boundary()) &&
               (p1.is_on_y_boundary() && p2.is_on_y_boundary())) ||
              ((p1.is_on_x_boundary() && p2.is_on_x_boundary()) && (y1 == y2)) ||
              ((x1 == x2)) && (p1.is_on_y_boundary() && p2.is_on_y_boundary()) ||
              ((x1 == x2) && (y1 == y2)));
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
      return (operator()(xcv1.left(), xcv2.left()) &&
              operator()(xcv1.right(), xcv2.right()))
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

    /*! Obtains the parameter space at a point along the x-axis.
     * \param p the point.
     * \return the parameter space at p.
     * \pre p does not lie on the x-boundary.
     */
    Arr_parameter_space operator()(const Point_2 p) const
    {
      CGAL_error();
      return ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_x_2 function object.
   * \return an object of type Parameter_space_in_x_2.
   */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(); }

  /*! A function object that obtains the parameter space of a geometric
   * entity along the y-axis
   */
  class Parameter_space_in_y_2 {
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
     *   ARR_TOP_BOUNDARY     - the arc reaches the y-boundary from below
     *                          at the curve ce end.
     * \pre xcv does not lie on the y-boundary.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const
    {
      CGAL_precondition(!m_traits.is_on_y_identification_2_object()(xcv));
      return (ce == ARR_MIN_END) ?
        (xcv.is_directed_top() ?
         xcv.source_parameter_space_in_x() : xcv.target_parameter_space_in_x()) :
        (xcv.is_directed_top() ?
         xcv.target_parameter_space_in_x() : xcv.source_parameter_space_in_x());
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     * \pre p does not lie on the y-boundary.
     */
    Arr_parameter_space operator()(const Point_2 p) const
    {
      CGAL_error();
      return ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_y_2 function object.
   * \return an object of type Parameter_space_in_y_2.
   */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }


  //! A functor that compares the x-coordinate of curve ends on the boundary of
  // the parameter space with curve ends and points.
  class Compare_x_on_boundary_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    //! The traits (in case it has state).
    const Traits& m_traits;

    //! Constructor
    // \param traits the traits (in case it has state)
    Compare_x_on_boundary_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    //! Compare the x-coordinate of a point with the x-coordinate of an
    // x-monotone curve end on the boundary.
    // \param point the point.
    // \param xcv the curve, the endpoint of which is compared.
    // \param ce the curve-end indicator -
    //            ARR_MIN_END - the minimal end of xc or
    //            ARR_MAX_END - the maximal end of xc.
    // \return the comparison result:
    //         SMALLER - x(p) < x(xcv, ce);
    //         EQUAL   - x(p) = x(xcv, ce);
    //         LARGER  - x(p) > x(xcv, ce).
    // \pre p does not lie on the x-boundary of the parameter space.
    // \pre the ce end of the curve xcv lies on the y-boundary of the
    //      parameter space.
    // \pre xcv does not lie on the x-boundary of the parameter space.
    Comparison_result operator()(const Point_2& point,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(point));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv));

      if (((ce == ARR_MIN_END) && xcv.is_directed_top()) ||
          ((ce == ARR_MAX_END) && !xcv.is_directed_top()))
      {
        CGAL_precondition(xcv.source_parameter_space_in_y ==
                          ARR_BOTTOM_BOUNDARY);
        Arr_parameter_space ps_x = xcv.source_parameter_space_in_x();
        return (ps_x == ARR_LEFT_BOUNDARY) ? LARGER :
          ((ps_x == ARR_RIGHT_BOUNDARY) ? SMALLER :
           (m_traits->compare_x_object()(p, xcv.source())));
      }
      // (ce == ARR_MIN_END) && !(xcv.is_directed_top() ||
      // (ce == ARR_MAX_END) && (xcv.is_directed_top() ||
      CGAL_precondition(xcv.target_parameter_space_in_y ==
                        ARR_BOTTOM_BOUNDARY);
      Arr_parameter_space ps_x = xcv.target_parameter_space_in_x();
      return (ps_x == ARR_LEFT_BOUNDARY) ? LARGER :
        ((ps_x == ARR_RIGHT_BOUNDARY) ? SMALLER :
         (m_traits->compare_x_object()(p, xcv.source())));
    }

    //! Compare the x-coordinates of 2 arc ends near the boundary of the
    // parameter space.
    // \param xcv1 the first arc.
    // \param ce1 the first arc end indicator -
    //            ARR_MIN_END - the minimal end of xcv1 or
    //            ARR_MAX_END - the maximal end of xcv1.
    // \param xcv2 the second arc.
    // \param ce2 the second arc end indicator -
    //            ARR_MIN_END - the minimal end of xcv2 or
    //            ARR_MAX_END - the maximal end of xcv2.
    // \return the second comparison result:
    //         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
    //         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
    //         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
    // \pre the ce1 end of the arc xcv1 lies on the y-boundary.
    // \pre the ce2 end of the arc xcv2 lies on the y-boundary.
    // \pre xcv1 does not lie in the x-identification curve.
    // \pre xcv2 does not lie in the x- identification curve.
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    {
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv1));
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(xcv2));

      Arr_parameter_space ps1_x;
      const Point_2* p1 = NULL;
      if (((ce1 == ARR_MIN_END) && xcv1.is_directed_top()) ||
          ((ce1 == ARR_MAX_END) && !xcv1.is_directed_top()))
      {
        Arr_parameter_space ps1_x = xcv1.source_parameter_space_in_x();
        p1 = &(xcv1.source());
      }
      else {
        Arr_parameter_space ps1_x = xcv1.target_parameter_space_in_x();
        p1 = &(xcv1.target());
      }

      Arr_parameter_space ps2_x;
      const Point_2* p1 = NULL;
      if (((ce2 == ARR_MIN_END) && xcv2.is_directed_top()) ||
            ((ce2 == ARR_MAX_END) && !xcv2.is_directed_top()))
      {
        Arr_parameter_space ps2_x = xcv1.source_parameter_space_in_x();
        p2 = &(xcv2.source());
      }
      else {
        Arr_parameter_space ps2_x = xcv1.target_parameter_space_in_x();
        p2 = &(xcv2.target());
      }

      if ((ps1_x == ARR_INTERIOR) && (ps2_x == ARR_INTERIOR))
        return m_traits->compare_x_object()(*p1, *p2);
      if (ps1_x == ps2_x) return EQUAL;
      return (ps1_x == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
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
      return m_traits->compare_x_object()(p1, p2);
    }
  };

  /*! Obtain a Compare_x_on_boundary_2 function object.
   * \return an object of type Compare_x_on_boundary_2.
   */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
  { return Compare_x_on_boundary_2(*this); }


  /*! A functor that compares the x-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_boundary_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
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
     * \pre the the x-coordinates of xcv1 and xcv2 at their ce end are equal.
     * \pre xcv1 does not lie on the x-boundary.
     * \pre xcv2 does not lie on the x-boundary.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    {
      //! \todo
      return EQUAL;
    }
  };

  /*! Obtain a Compare_x_near_boundary_2 function object.
   * \return an object of type Compare_x_near_boundary_2.
   */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  { return Compare_x_near_boundary_2(*this); }


  /*! A functor that compares the y-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_near_boundary_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Compare the y-coordinates of 2 x-monotone curves at their ends near
     * the boundary of the parameter space.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the arc end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the arcs xcv1 and xcv2 lie either on the left
     *      boundary or on the right boundary of the parameter space (implying
     *      that they cannot be vertical).
     * There is no horizontal identification curve!
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    {
      //! \todo
      CGAL_precondition(!xcv1.is_degenerate());
      CGAL_precondition(!xcv2.is_degenerate());
      return EQUAL;
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

    /*! The traits (in case it has state) */
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
      return m_traits->compare_y_object()(p1, p2);
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
    bool operator()(const Point_2& p) const { return p._is_on_x_boundary(); }

    /*! Determine whether an x-monotone curve lies on the y-boundary.
     * \param xcv the curve.
     * \return a Boolean indicating whether xcv lies on the y-boundary.
     */
    bool operator()(const X_monotone_curve_2& xcv) const
    { return (operator()(xcv.source()) && operator()(xcv.target())); }
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

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
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
    { return (operator()(xcv.source()) && operator()(xcv.target())); }
  };

  /*! Obtain a Is_on_y_identification_2 function object.
   * \return an object of type Is_on_y_identification_2.
   */
  Is_on_y_identification_2 is_on_y_identification_2_object() const
  { return Is_on_y_identification_2(*this); }
  //@}

  /// \name Functor definitions for supporting intersections.
  //@{

  /*! A functor that divides an arc into x-monotone arcs. That are, arcs that
   * do not cross the identification arc.
   */
  class Make_x_monotone_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Make_x_monotone_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Cut the given curve into x-monotone subcurves and insert them into the
     * given output iterator. As spherical_arcs are always x_monotone, only one
     * object will be contained in the iterator.
     * \param xc the curve.
     * \param oi the output iterator, whose value-type is Object. The output
     *           object is a wrapper of either an X_monotone_curve_2, or - in
     *           case the input spherical_arc is degenerate - a Point_2 object.
     * \return the past-the-end iterator.
     */
    template<typename OutputIterator>
    OutputIterator operator()(const Curve_2& c, OutputIterator oi) const
    {
      //!\todo
      if (c.is_degenerate()) {
        // The spherical_arc is a degenerate point - wrap it with an object:
        *oi++ = make_object(c.right());
        return oi;
      }

      if (c.is_x_monotone()) {
        // The spherical arc is monotone - wrap it with an object:
        // *oi++ = make_object(X_monotone_curve_2(c));
        const X_monotone_curve_2* xc = &c;
        *oi++ = make_object(*xc);
        return oi;
      }

      // if (c.is_full()) {
      // }

      //! \todo
      return oi;
    }
  };

  /*! Obtain a Make_x_monotone_2 function object.
   * \return an object of type Make_x_monotone_2.
   */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(*this); }

  /*! A functor that splits an x-monotone arc at a directional point. */
  class Split_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Split_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Split a given x-monotone curve at a given point into two sub-curves.
     * \param xc the curve to split
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
      //! \todo
      CGAL_precondition(!xcv.is_degenerate());
    }
  };

  /*! Obtain a Split_2 function object. */
  Split_2 split_2_object() const { return Split_2(*this); }

  /*! A functor that computes intersections between x-monotone arcs. */
  class Intersect_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Intersect_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Find the intersections of the two given curves and insert them into the
     * given output iterator. As two spherical_arcs may itersect only once,
     * only a single intersection will be contained in the iterator.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param oi the output iterator.
     * \return the past-the-end iterator.
     * \pre xcv1 and xcv2 are not degenerate
     */
    template<typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv1,
                              const X_monotone_curve_2& xcv2,
                              OutputIterator oi) const
    {
      CGAL_precondition(!xcv1.is_degenerate());
      CGAL_precondition(!xcv2.is_degenerate());

      //! \todo
      return oi;
    }
  };

  /*! Obtain an Intersect_2 function object. */
  Intersect_2 intersect_2_object() const { return Intersect_2(*this); }

  /*! A functor that tests whether two x-monotone arcs can be merged. */
  class Are_mergeable_2 {
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Are_mergeable_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \return true if the two curves are mergeable; false otherwise.
     * Two arcs are mergeable if:
     * 1. they are supported by the same plane, and
     * 2. share a common endpoint that is not on the identification arc
     */
    bool operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2) const
    {
      if (xcv1.is_empty() || xcv2.is_empty()) return true;
      if ((xcv1.is_x_full() || xcv1.is_meridian()) &&
          (xcv2.is_x_full() || xcv2.is_meridian())) return false;

      //! \todo
      return false;
    }
  };

  /*! Obtain an Are_mergeable_2 function object.
   * \return an object of type Are_mergeable_2.
   */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(*this); }

  /*! A functor that merges two x-monotone arcs into one */
  class Merge_2 {
  protected:
    typedef Arr_flat_torus_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Merge two given x-monotone curves into a single curve (spherical_arc).
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param xcv Output: the merged curve.
     * \pre the two curves are mergeable.
     */
    void operator()(const X_monotone_curve_2& xcv1,
                    const X_monotone_curve_2& xcv2,
                    X_monotone_curve_2& xcv) const
    {
      CGAL_precondition(m_traits->are_mergeable_2_object()(xcv1, xcv2) == true);

      if (xcv1.is_degenerate() || xcv1.is_empty()) {
        xc = xcv2;
        return;
      }

      if (xcv2.is_degenerate() || xcv2.is_empty()) {
        xc = xcv1;
        return;
      }
      //! \todo
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
    /*! Return an approximation of a point coordinate.
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

    //! The traits.
    const Traits& m_traits;

    //! Constructor
    // \param traits the traits (in case it has state)
    onstruct_x_monotone_curve_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_flat_torus_traits_2<Kernel>;

  public:
    /*! Return an x-monotone curve connecting the two given endpoints.
     * \param p1 the first point.
     * \param p2 the second point.
     * \pre p1 and p2 must not be the same.
     * \return an x-monotone torusial arc connecting p1 and p2.
     */
    X_monotone_curve_2 operator()(const Point_2& source, const Point_2& target,
                                  bool is_vertical,
                                  bool is_directed_right,
                                  bool is_directed_top,
                                  bool is_degenerate = false,
                                  bool is_empty = false,
                                  bool m_is_x_full = false) const
    {
      return X_monotone_curve_2(source, target, is_vertical,
                                is_directed_right, is_directed_top,
                                is_x_full, is_degenerate, is_empty);
    }

    //! Return an x-monotone curve connecting the two given endpoints.
    // \param p1 the first point.
    // \param p2 the second point.
    // \pre p1 and p2 must not be the same.
    // \return an x-monotone torusial arc connecting p1 and p2.
    X_monotone_curve_2 operator()(const Point_2& source,
                                  const Point_2& target) const
    {
      CGAL_precondition(!kernel.equal_2_object()(p1, p2));
      bool is_degenerate(false);
      bool is_empty(false);
      bool m_is_full(((source.x() == 0) && (target.x() == 1)) ||
                     ((target.x() == 1) && (source.x() == 0)));
      bool is_directed_right((source.x() < target.x()) ||
                             ((source.x() == target.x()) &&
                              (source.y() < target.y())));
      bool is_directed_top((source.y() < target.y()) ||
                             ((source.y() == target.y()) &&
                              (source.x() < target.x())));
      return X_monotone_curve_2(source, target, is_vertical,
                                is_directed_right, is_directed_top,
                                is_x_full, is_degenerate, is_empty);
    }
  };

  /*! Obtain a Construct_x_monotone_curve_2 function object.
   * \return an object of type Construct_x_monotone_curve_2.
   */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(*this); }
  //@}

  /// \name Functor definitions for the Boolean set-operation traits.
  //@{
  class Compare_endpoints_xy_2 {
  public:
    /*! Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param xc the curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv)
    { return (xc.is_directed_right()) ? SMALLER : LARGER; }
  };

  /*! Obtain a Compare_endpoints_xy_2 function object.
   * \return an object of type Compare_endpoints_xy_2.
   */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }

  class Construct_opposite_2 {
  public:
    /*! Construct an opposite x-monotone (with swapped source and target).
     * \param xc the curve.
     * \return the opposite curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xc)
    { return xc.opposite(); }
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
      CGAL_precondition(!m_traits.is_on_x_identification_2_object()(p));

      const FT& px = p.x();
      const FT& left_x = xcv.is_directed_right() ? xcv.source_x() : target_x();
      const FT& right_x = xcv.is_directed_right() ? xcv.target_x() : source_x();

      return ((left_x <= px) && (px <= right_x));
    }
  };

  /*! Obtain an Is_in_x_range function object.
   * \return an object of type Is_in_x_range.
   */
  Is_in_x_range is_in_x_range_object() const { return Is_in_x_range(); }
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
    Point_2(x, y),
    m_on_x_boundary((x == 0) || (x == 1))
    m_on_y_boundary((y == 0) || (y == 1))
  {}

  /*! Constructor from a point.
   * \param x the x-coordinate.
   * \param y the y-coordinate.
   */
  Arr_point_on_flat_torus_3(const Point_2& p) :
    Point_2(p),
    m_on_x_boundary((p.x() == 0) || (p.x() == 1))
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
  { return is_on_ve_bountical_boundary() || is_on_y_boundary() }
};

//! A Representation of an x-monotone great circular arc embedded on a torus,
// as used by the Arr_flat_torus_traits_2 traits-class template.
// An x-monotone great circular arc cannot cross boundary of the parameter
// space.
// \todo At this point such an arc cannot have an angle of 180 degrees.
template <typename Kernel_>
class Arr_x_monotone_curve_on_flat_torus_3 {
public:
  typedef Kernel_                                   Kernel;

  typedef typename Kernel::FT                       FT;

protected:
  // For some reason compilation under Windows fails without the qualifier
  typedef CGAL::Arr_point_on_flat_torus_3<Kernel>   Arr_point_on_flat_torus_3;

  //! The source point of the arc.
  Arr_point_on_flat_torus_3 m_source;

  //! The target point of the arc.
  Arr_point_on_flat_torus_3 m_target;

  //! The arc is vertical.
  bool m_is_vertical;

  //! Target is xy-lexicographically larger than source.
  bool m_is_directed_right;

  //! Target is yx-lexicographically larger than source.
  bool m_is_directed_top;

  //! The arc is a full circle.
  bool m_is_full;

  // The arc is degenerate---it consists of a single point.
  bool m_is_degenerate;

  //! The arc is empty.
  bool m_is_empty;

public:
  //! Default constructor - constructs an empty arc.
  Arr_x_monotone_curve_on_flat_torus_3() :
    m_is_vertical(false),
    m_is_directed_right(false),
    m_is_directed_top(false),
    m_is_full(false),
    m_is_degenerate(false),
    m_is_empty(true)
  {}

  //! Constructor
  // \param src the source point of the arc
  // \param trg the target point of the arc
  // \param num_horizontal_revolutions the number of horizontal of the arc
  // \param num_vertical_revolutions the number of vertical of the arc
  // \param is_vertical is the arc vertical ?
  // \param is_directed_right is the curve directed from left to right?
  // \param is_directed_top is the curve directed from bottom to top?
  // \param is_full is the arc a full circle?
  // \param is_degenerate is the arc degenerate (single point)?
  Arr_x_monotone_curve_on_flat_torus_3(const Arr_point_on_flat_torus_3& src,
                                       const Arr_point_on_flat_torus_3& trg,
                                       bool is_vertical,
                                       bool is_directed_right,
                                       bool is_directed_top,
                                       bool is_full = false,
                                       bool is_degenerate = false,
                                       bool is_empty = false) :
    m_source(src),
    m_target(trg),
    m_is_vertical(is_vertical),
    m_is_directed_right(is_directed_right),
    m_is_directed_top(is_directed_top),
    m_is_full(is_full),
    m_is_degenerate(is_degenerate),
    m_is_empty(is_empty)
  {}

  /*! Set the source endpoint direction.
   * \param p the endpoint to set.
   */
  void set_source(const Arr_point_on_flat_torus_3& p) { m_source = p; }

  /*! Set the target endpoint direction.
   * \param p the endpoint to set.
   */
  void set_target(const Arr_point_on_flat_torus_3& p) { m_target = p; }

  void set_is_vertical(bool flag) { m_is_vertical = flag; }
  void set_is_directed_right(bool flag) { m_is_directed_right = flag; }
  void set_is_directed_top(bool flag) { m_is_directed_top = flag; }
  void set_is_full(bool flag) { m_is_full = flag; }
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

  /*! Determine whether the curve is vertical */
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
  bool is_x_full() const { return m_is_x_full; }

  /*! Determine whether the curve is degenerate */
  bool is_degenerate() const { return m_is_degenerate; }

  /*! Determine whether the curve is degenerate */
  bool is_empty() const { return m_is_empty; }

  /*! Obtain the parameter space in x of the source point.
   */
  Arr_parameter_space source_parameter_space_in_x() const
  {
    return (!source.is_on_x_boundary()) ? ARR_INTERIOR :
      (is_directed_right() ? ARR_LEFT_BOUNDARY : ARR_RIGHT_BOUNDARY) ;
  }

  /*! Obtain the parameter space in y of the source point.
   */
  Arr_parameter_space source_parameter_space_in_y() const
  {
    return (!source.is_on_y_boundary()) ? ARR_INTERIOR :
      (is_directed_top() ? ARR_BOTTOM_BOUNDARY : ARR_TOP_BOUNDARY) ;
  }

  /*! Obtain the parameter space in x of the target point.
   */
  Arr_parameter_space target_parameter_space_in_x() const
  {
    return (!target.is_on_x_boundary()) ? ARR_INTERIOR :
      (is_directed_right() ? ARR_RIGHT_BOUNDARY : ARR_LEFT_BOUNDARY) ;
  }

  /*! Obtain the parameter space in y of the target point.
   */
  Arr_parameter_space target_parameter_space_in_y() const
  {
    return (!target.is_on_y_boundary()) ? ARR_INTERIOR :
      (is_directed_top() ? ARR_TOP_BOUNDARY : ARR_BOTTOM_BOUNDARY) ;
  }

  /*! Obtain the x-coordinate of the source point.
   */
  const FT& source_x() const
  {
    return (!source.is_on_x_boundary()) ? source.x() :
      (is_directed_right() ? 0 : 1) ;
  }

  /*! Obtain the y-coordinate of the source point.
   */
  const FT& source_y() const
  {
    return (!source.is_on_y_boundary()) ? source.y() :
      (is_directed_top() ? 0 : 1) ;
  }

  /*! Obtain the x-coordinate of the target point.
   */
  const FT& target_x() const
  {
    return (!target.is_on_x_boundary()) ? target.x() :
      (is_directed_right() ? 1 : 0) ;
  }

  /*! Obtain the y-coordinate of the target point.
   */
  const FT& target_y() const
  {
    return (!target.is_on_y_boundary()) ? target.y() :
      (is_directed_top() ? 1 : 0) ;
  }

  /*! Determine whether the curve lie on the boundary. */
  bool is_on_boundary() const
  {
    return (source_parameter_space_in_x == target_parameter_space_in_x) ||
      (source_parameter_space_in_y == target_parameter_space_in_y);
  }

  /*! Flip the arc (swap it source and target) */
  Arr_x_monotone_curve_on_flat_torus_3 opposite() const
  {
    Arr_x_monotone_curve_on_flat_torus_3 opp;
    opp.m_source = this->m_target;
    opp.m_target = this->m_source;
    opp.m_is_directed_right = !(this->is_directed_right());
    opp.m_is_directed_top = !(this->is_directed_top());
    opp.m_is_vertical = this->is_vertical();
    opp.m_is_full = this->is_full();
    opp.m_is_degenerate = this->is_degenerate();
    opp.m_is_empty = this->is_empty();
    return opp;
  }
};

//! A representation of a geodesic arc embedded on a flat torus,
// used by the Arr_flat_torus_traits_2 traits-class template.
// An arc is uniqely represented by two endpoints, the source s and the
// target t. The points of the arc are the locus of points visited when
// moving from the source s toward the target t in a straight line in the
// parameter space.
template <typename Kernel_>
class Arr_curve_on_flat_torus_3 :
  public Arr_x_monotone_curve_on_flat_torus_3<Kernel_> {
public:
  typedef Kernel_                                      Kernel;

protected:
  typedef Arr_x_monotone_curve_on_flat_torus_3<Kernel> Base;

  typedef typename Base::Arr_point_on_flat_torus_3     Arr_point_on_flat_torus_3;

private:
  //! The number of revolutions in the horizontal direction.
  size_t m_num_horizontal_revolutions;

  //! The number of revolutions in the vertical direction.
  size_t m_num_vertical_revolutions;

public:
  //! Default constructor - constructs an empty arc.
  Arr_curve_on_flat_torus_3() :
    Base(),
    m_num_horizontal_revolutions(0),
    m_num_vertical_revolutions(0)
  {}

  /*! Constructor
   * \param src the source point of the arc
   * \param trg the target point of the arc
   * \param xxx
   * \param is_x_monotone is arc  x-monotone ?
   * \param is_vertical is the arc vertical ?
   * \param is_directed_right is the arc directed from left to right?
   * \param is_full is the arc a full (great) circle?
   * \param is_degenerate is the arc degenerate (single point)?
   * \pre plane contains the origin
   * \pre plane contains src
   * \pre plane contains trg
   */
  Arr_curve_on_flat_torus_3(const Arr_extended_direction_3& src,
                            const Arr_extended_direction_3& trg,
                            int num_horizontal_revolutions,
                            int num_vertical_revolutions,
                            bool is_vertical,
                            bool is_directed_right,
                            bool is_full = false,
                            bool is_degenerate = false,
                            bool is_empty = false) :
    Base(src, trg, is_vertical, is_directed_right, is_directed_top,
         is_full, is_degenerate, is_empty),
    m_num_horizontal_revolutions(num_horizontal_revolutions),
    m_num_vertical_revolutions(num_vertical_revolutions)
  {}

  //! Determines the number of horizontal_revolutions of the arc.
  // \return the number of horizontal_revolutions.
  int number_of_horizontal_revolutions() const
  { return m_num_horizontal_revolutions; }

  //! Determines the number of vertical_revolutions of the arc.
  // \return the number of vertical_revolutions.
  int number_of_vertical_revolutions() const
  { return m_num_vertical_revolutions; }
};

/*! Inserter for the flat-torus point used by the traits-class.
 */
template <typename Kernel_, typename InputStream>
InputStream& operator>>(InputStream& os,
                        const Arr_point_on_flat_torus_3<Kernel_>& point)
{
  const Kernel::Point_2* p = &point;
  os << *p
  return os;
}

/*! Inserter for the flat-torus x-monotone curve used by the traits-class
 */
template <typename Kernel, typename OutputStream>
OutputStream&
operator<<(OutputStream& os,
           const Arr_x_monotone_curve_on_flat_torus_3<Kernel_>& xcv)
{
  os << "(";
  os << xcv.source() << ", " << xcv.target();
#if defined(CGAL_ARR_CURVE_ON_FLAT_TORUS_DETAILS)
#endif
  os << ")";
  return os;
}

/*! Inserter for the flat-torus curve used by the traits-class */
template <typename Kernel, typename OutputStream>
OutputStream&
operator<<(OutputStream& os, const Arr_curve_on_flat_torus_3<Kernel>& cv)
{
  os << "(";
  os << xcv.source() << ", " << xcv.target();
#if defined(CGAL_ARR_CURVE_ON_FLAT_TORUS_DETAILS)
#endif
  os << ")";
  return os;
}

/*! Extractor for the flat-torus point used by the traits-class */
template <typename Kernel_, typename InputStream>
InputStream& operator>>(InputStream& is,
                        Arr_point_on_flat_torus_3<Kernel_>& point)
{
  Kernel::Point_2 p;
  is >> p;
  point(p);
  return is;
}

/*! Extractor for the flat-torus x-monotone curve used by the traits-class */
template <typename Kernel_, typename InputStream>
InputStream&
operator>>(InputStream& is, Arr_x_monotone_curve_on_flat_torus_3<Kernel_>& xcv)
{
  Point_2 source, target;
  is >> source;
  is >> target;
  xcv(source, target);
  return is;
}

} //namespace CGAL

#endif
