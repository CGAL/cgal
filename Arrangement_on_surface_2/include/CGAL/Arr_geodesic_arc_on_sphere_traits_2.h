// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_GEODESIC_ARC_ON_SPHERE_TRAITS_2_H
#define CGAL_ARR_GEODESIC_ARC_ON_SPHERE_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

// #define CGAL_FULL_X_MONOTONE_GEODESIC_ARC_ON_SPHERE_IS_SUPPORTED 1

/*! \file
 * A class that handles great circular arcs embedded on spheres suitable
 * as a geometry traits class for the arrangement on surface package.
 */

#include <fstream>

#include <boost/variant.hpp>

#include <CGAL/config.h>
#include <CGAL/tags.h>
#include <CGAL/tss.h>
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

template <typename Kernel> class Arr_x_monotone_geodesic_arc_on_sphere_3;
template <typename Kernel> class Arr_geodesic_arc_on_sphere_3;
template <typename Kernel> class Arr_extended_direction_3;

/*! A traits class-template for constructing and maintaining arcs of great
 * circles embedded on spheres. It is parameterized from a (linear) geometry
 * kernel, which it also derives from
 */
template <typename T_Kernel>
class Arr_geodesic_arc_on_sphere_traits_2 : public T_Kernel {
  friend class Arr_x_monotone_geodesic_arc_on_sphere_3<T_Kernel>;
  friend class Arr_geodesic_arc_on_sphere_3<T_Kernel>;
  friend class Arr_extended_direction_3<T_Kernel>;

public:
  typedef T_Kernel                              Kernel;

  // Category tags:
  typedef Tag_true                              Has_left_category;
  typedef Tag_true                              Has_merge_category;
  typedef Tag_false                             Has_do_intersect_category;

  typedef Arr_identified_side_tag               Left_side_category;
  typedef Arr_contracted_side_tag               Bottom_side_category;
  typedef Arr_contracted_side_tag               Top_side_category;
  typedef Arr_identified_side_tag               Right_side_category;

  /*! Default constructor */
  Arr_geodesic_arc_on_sphere_traits_2(){}

protected:
  typedef typename Kernel::FT                   FT;

  typedef typename Kernel::Direction_3          Direction_3;
  typedef typename Kernel::Vector_3             Vector_3;
  typedef typename Kernel::Direction_2          Direction_2;
  typedef typename Kernel::Vector_2             Vector_2;

  /*! Obtain the possitive (north) pole
   * \return the possitive (north) pole
   */
  inline static const Direction_3& pos_pole()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_3(Direction_3, d, 0, 0, 1);
    return d;
  }

  /*! Obtain the negative (south) pole
   * \return the negative (south) pole
   */
  inline static const Direction_3& neg_pole()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_3(Direction_3, d, 0, 0, -1);
    return d;
  }

  /*! Obtain the intersection of the identification arc and the xy plane.
   * By default, it is the vector directed along the negative x axis
   * (x = -infinity).
   * \return the intersection of the identification arc and the xy plane.
   */
  inline static const Direction_2& identification_xy()
  {
#if (CGAL_IDENTIFICATION_XY == CGAL_X_MINUS_1_Y_0)
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(Direction_2, d, -1, 0);

#elif (CGAL_IDENTIFICATION_XY == CGAL_X_MINUS_8_Y_6)
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(Direction_2, d, -8, 6);

#elif (CGAL_IDENTIFICATION_XY == CGAL_X_MINUS_11_Y_7)
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(Direction_2, d, -11, 7);
#else
#error CGAL_IDENTIFICATION_XY is not defined
#endif
    return d;
  }

  /*! Obtain the normal of the plane that contains the identification arc.
   * By default, it is the vector directed along the positive y axis
   * (y = infinity).
   * \return the normal of the plane that contains the identification arc.
   */
  inline static const Direction_3& identification_normal()
  {
#if (CGAL_IDENTIFICATION_XY == CGAL_X_MINUS_1_Y_0)
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_3(Direction_3, d, 0, 1, 0);
#elif (CGAL_IDENTIFICATION_XY == CGAL_X_MINUS_8_Y_6)
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_3(Direction_3, d, 6, 8, 0);
#elif (CGAL_IDENTIFICATION_XY == CGAL_X_MINUS_11_Y_7)
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_3(Direction_3, d,7, 11, 0);
#else
#error CGAL_IDENTIFICATION_XY is not defined
#endif
    return d;
  }

  /*! Obtain the 2D direction directed along the negative x axis
   * \return the direction directed at x = -infinity
   */
  inline static const Direction_2& neg_x_2()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(Direction_2, d, -1, 0);
    return d;
  }

  /*! Obtain the 2D direction directed along the negative y axis
   * \return the direction directed at y = -infinity
   */
  inline static const Direction_2& neg_y_2()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(Direction_2, d, 0, -1);
    return d;
  }

  /*! Obtain the sign of the x-coordinate of a direction in space
   * \param d the direction in space
   * \return the sign of the x-coordinate of d
   */
  inline static Sign x_sign(Direction_3 d) { return CGAL::sign(d.dx()); }

  /*! Obtain the sign of the y-coordinate of a direction in space
   * \param d the direction in space
   * \return the sign of the y-coordinate of d
   */
  inline static Sign y_sign(Direction_3 d) { return CGAL::sign(d.dy()); }

  /*! Obtain the sign of the z-coordinate of a direction in space
   * \param d the direction in space
   * \return the sign of the z-coordinate of d
   */
  inline static Sign z_sign(Direction_3 d) { return CGAL::sign(d.dz()); }

  typedef Direction_2 (*Project)(const Direction_3& d) ;

  /*! Project a 3D direction onto the xy-plane
   * \param d the 3D direction
   * \return the projection onto the xy-plane
   */
  inline static Direction_2 project_xy(const Direction_3& d)
  { return Direction_2(d.dx(), d.dy()); }

  /*! Project a 3D direction onto the yz-plane
   * \param d the 3D direction
   * \return the projection onto the yz-plane
   */
  inline static Direction_2 project_yz(const Direction_3& d)
  { return Direction_2(d.dy(), d.dz()); }

  /*! Project a 3D direction onto the zx-plane
   * \param d the 3D direction
   * \return the projection onto the xz-plane
   */
  inline static Direction_2 project_xz(const Direction_3& d)
  { return Direction_2(d.dx(), d.dz()); }

  /*! Project a 3D direction onto the yz-plane and reflect the y coordinate
   * \param d the 3D direction
   * \return the projection onto the yz-plane
   */
  inline static Direction_2 project_minus_yz(const Direction_3& d)
  { return Direction_2(-d.dy(), d.dz()); }

  /*! Project a 3D direction onto the zx-plane and reflect the x coordinate
   * \param d the 3D direction
   * \return the projection onto the xz-plane
   */
  inline static Direction_2 project_minus_xz(const Direction_3& d)
  { return Direction_2(-d.dx(), d.dz()); }

  /*! Compare the relative position of a direction and a plane given by its
   * normal.
   * \param normal the direction of the plane.
   * \param dir the direction.
   */
  inline Oriented_side oriented_side(const Direction_3& normal,
                                     const Direction_3 dir) const
  {
    FT dot = normal.vector() * dir.vector();
    return CGAL::sign(dot);
  }

  /*! Compute the orientation of two directions.
   * \param d1 the first direction.
   * \param d2 the second direction.
   * \return the relative orientation of d1 and d2.
   */
  static inline Orientation orientation(const Direction_2& d1,
                                        const Direction_2& d2)
  {
    Kernel kernel;
    return kernel.orientation_2_object()(d1.vector(), d2.vector());
  }

  /*! Determined whether a direction is contained in a plane given by its
   * normal
   * \param normal the direction of the 3D plane.
   * \param dir the 3D direction.
   * \return true if dir is contained in plane; false otherwise.
   * \pre the plane contains the origin.
   */
  inline bool has_on(const Direction_3& normal, const Direction_3& dir) const
  {
    FT dot = normal.vector() * dir.vector();
    return CGAL::sign(dot) == ZERO;
  }

public:
  /*! Compare two endpoint directions by v.
   * \param d1 the first enpoint direction.
   * \param d2 the second endpoint direction.
   * \return SMALLER - v(d1) < v(d2);
   *         EQUAL   - v(d1) = v(d2);
   *         LARGER  - v(d1) > v(d2).
   */
  inline Comparison_result compare_y(const Direction_3& d1,
                                     const Direction_3& d2) const
  {
    Vector_3 v1 = d1.vector();
    Vector_3 v2 = d2.vector();

    FT dot_p1 = v1.z();
    FT dot_p2 = v2.z();

    Sign s1 = CGAL::sign(dot_p1);
    Sign s2 = CGAL::sign(dot_p2);

    if (s1 != s2) {
      if (s1 == NEGATIVE) return SMALLER;
      if (s1 == POSITIVE) return LARGER;
      if (s2 == NEGATIVE) return LARGER;
      if (s2 == POSITIVE) return SMALLER;
    }
    if (s1 == ZERO) return EQUAL;

    FT norm1 = v1 * v1;
    FT norm2 = v2 * v2;

    return (s1 == POSITIVE) ?
      CGAL::compare(dot_p1 * dot_p1 * norm2, dot_p2 * dot_p2 * norm1) :
      CGAL::compare(dot_p2 * dot_p2 * norm1, dot_p1 * dot_p1 * norm2);
  }

  /*! Compare two directions contained in the xy plane by u.
   * \param d1 the first direction.
   * \param d2 the second direction.
   * \return SMALLER - u(d1) < u(d2);
   *         EQUAL   - u(d1) = u(d2);
   *         LARGER  - u(d1) > u(d2).
   */
  inline Comparison_result compare_x(const Direction_2& d1,
                                     const Direction_2& d2) const
  {
    const Kernel* kernel = this;
    if (kernel->equal_2_object()(d1, d2)) return EQUAL;
    const Direction_2& d = identification_xy();
    return (kernel->counterclockwise_in_between_2_object()(d, d1, d2)) ?
      LARGER : SMALLER;
  }

  /*! Compare two endpoint directions by u.
   * \param d1 the first enpoint direction.
   * \param d2 the second endpoint direction.
   * \return SMALLER - u(d1) < u(d2);
   *         EQUAL   - u(d1) = u(d2);
   *         LARGER  - u(d1) > u(d2).
   * \pre d1 does not coincide with any pole.
   * \pre d2 does not coincide with any pole.
   */
  inline Comparison_result compare_x(const Direction_3& d1,
                                     const Direction_3& d2) const
  {
    // Compare the projections onto the xy plane:
    Direction_2 d1_2 = project_xy(d1);
    Direction_2 d2_2 = project_xy(d2);
    return compare_x(d1_2, d2_2);
  }

  /*! Compare two endpoint directions lexigoraphically: by u, then by v.
   * \param d1 the first enpoint direction.
   * \param d2 the second endpoint direction.
   * \return SMALLER - u(d1) < u(d2);
   *         SMALLER - u(d1) = u(d2) and v(d1) < v(d2);
   *         EQUAL   - u(d1) = u(d2) and v(d1) = v(d2);
   *         LARGER  - u(d1) = u(d2) and v(d1) > v(d2);
   *         LARGER  - u(d1) > u(d2).
   * \pre d1 does not lie on the discontinuity arc.
   * \pre d2 does not lie on the discontinuity arc.
   */
  inline Comparison_result compare_xy(const Direction_3& d1,
                                      const Direction_3& d2) const
  {
    Comparison_result res = compare_x(d1, d2);
    if (res == EQUAL) return compare_y(d1, d2);
    return res;
  }

public:
  // Traits objects
  typedef Arr_extended_direction_3<Kernel>                Point_2;
  typedef Arr_x_monotone_geodesic_arc_on_sphere_3<Kernel> X_monotone_curve_2;
  typedef Arr_geodesic_arc_on_sphere_3<Kernel>            Curve_2;
  typedef unsigned int                                    Multiplicity;

public:
  /// \name Basic functor definitions
  //@{

  /*! A functor that compares the x-coordinates of two directional points */
  class Compare_x_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_x_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Compare the x-coordinates of two directional points.
     * \param p1 the first directional point.
     * \param p2 the second directional point.
     * \return SMALLER - x(p1) < x(p2);
     *         EQUAL   - x(p1) = x(p2);
     *         LARGER  - x(p1) > x(p2).
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
    */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      CGAL_precondition(p1.is_no_boundary());
      CGAL_precondition(p2.is_no_boundary());

      return m_traits->compare_x(p1, p2);
    }
  };

  /*! Obtain a Compare_x_2 function object */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(this); }

  /*! A functor that compares two directional points lexigoraphically:
   * by x, then by y.
   */
  class Compare_xy_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_xy_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Compare two directional points lexigoraphically: by x, then by y.
     * \param p1 the first enpoint directional point.
     * \param p2 the second endpoint directional point.
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
      CGAL_precondition(p1.is_no_boundary());
      CGAL_precondition(p2.is_no_boundary());

      return m_traits->compare_xy(p1, p2);
    }
  };

  /*! Obtain a Compare_xy_2 function object */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(this); }

  /*! A functor that obtain the left endpoint of an x-monotone arc */
  class Construct_min_vertex_2 {
  public:
    /*! Obtain the left endpoint of and arc.
     * \param xc the arc.
     * \return the left endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& xc) const
    { return xc.left(); }
  };

  /*! Obtain a Construct_min_vertex_2 function object */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  /*! A functor that obtain the right endpoint of an x-monotone arc */
  class Construct_max_vertex_2 {
  public:
    /*! Obtain the right endpoint of an arc.
     * \param xc the arc.
     * \return the right endpoint.
     */
    const Point_2& operator()(const X_monotone_curve_2& xc) const
    { return xc.right(); }
  };

  /*! Obtain a Construct_max_vertex_2 function object */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  /*! A functor that checks whether an x-monotone arc is a vertical */
  class Is_vertical_2 {
  public:
    /*! Check whether a given arc is vertical.
     * \param xc the arc.
     * \return true if the curve is a vertical spherical_arc; false otherwise.
     * \pre the arc is not degenerate (consists of a single point)
     */
    bool operator()(const X_monotone_curve_2& xc) const
    {
      CGAL_precondition(!xc.is_degenerate());
      return xc.is_vertical();
    }
  };

  /*! Obtain an Is_vertical_2 function object */
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(); }

  /*! A functor that compares the y-coordinates of a directional point and
   * an x-monotone arc at the point x-coordinate
   */
  class Compare_y_at_x_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Return the location of the given point with respect to the input arc.
     * \param xc the arc.
     * \param p the point.
     * \return SMALLER - y(p) < xc(x(p)), i.e. the point is below the curve;
     *         EQUAL   - p lies on the curve.
     *         LARGER  - y(p) > xc(x(p)), i.e. the point is above the curve;
     * \pre p is not a construction point.
     * \pre p is in the x-range of xc.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xc) const
    {
      CGAL_precondition(!p.is_min_boundary() && !p.is_max_boundary());
      CGAL_precondition(xc.is_in_x_range(p));

      if (xc.is_vertical()) {
        // Compare the point with the left endpoint. If smaller, return SMALLER.
        // Otherwise, if EQUAL, return EQUAL.
        // Otherwise, compare with the right endpoint. If larger, return LARGER.
        // Otherwise, return EQUAL:
        if (!xc.left().is_min_boundary()) {
          Comparison_result cr = m_traits->compare_y(p, xc.left());
          if (cr != LARGER) return cr;
        }
        if (xc.right().is_max_boundary()) return EQUAL;
        Comparison_result cr = m_traits->compare_y(p, xc.right());
        return (cr == LARGER) ? LARGER : EQUAL;
      }

      // Compare the point to the underlying plane of xc:
      Oriented_side os = m_traits->oriented_side(xc.normal(), p);
      return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
        (xc.is_directed_right()) ?
        ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER) :
        ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER);
    }
  };

  /*! Obtain a Compare_y_at_x_2 function object */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(this); }

  /*! A functor that compares the y-coordinates of two x-monotone arcs
   * immediately to the left of their intersection directional point.
   */
  class Compare_y_at_x_left_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_left_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Compare the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \param p the intersection point.
     * \return the relative position of xc1 with respect to xc2 immdiately to
     *         the left of p: SMALLER, EQUAL, or LARGER.
     * \pre the point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \pre the arcs are not degenerate
     */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 const Point_2&
                                 CGAL_precondition_code(p)) const
    {
      CGAL_precondition(!xc1.is_degenerate());
      CGAL_precondition(!xc2.is_degenerate());
      CGAL_precondition(p == xc1.right());
      CGAL_precondition(p == xc2.right());

      // If Both arcs are vertical, they overlap:
      if (xc1.is_vertical() && xc2.is_vertical()) return EQUAL;
      if (xc1.is_vertical()) return SMALLER;
      if (xc2.is_vertical()) return LARGER;

      // Non of the arc is verticel. Thus, non of the endpoints coincide with
      // a pole.
      // Compare the y-coord. at the x-coord of the most right left-endpoint.
      const Point_2& l1 = xc1.left();
      const Point_2& l2 = xc2.left();
      if (!l1.is_no_boundary()) {
        // use l2 and xc1:
        Oriented_side os = m_traits->oriented_side(xc1.normal(), l2);
        return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
          (xc1.is_directed_right()) ?
          ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER) :
          ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER);
      }
      if (!l2.is_no_boundary()) {
        // use l1 and xc2:
        Oriented_side os = m_traits->oriented_side(xc2.normal(), l1);
        return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
          (xc2.is_directed_right()) ?
          ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER) :
          ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER);
      }

      if (m_traits->compare_xy(l1, l2) == SMALLER) {
        // use l2 and xc1:
        Oriented_side os = m_traits->oriented_side(xc1.normal(), l2);
        return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
          (xc1.is_directed_right()) ?
          ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER) :
          ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER);
      }
      // use l1 and xc2:
      Oriented_side os = m_traits->oriented_side(xc2.normal(), l1);
      return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
        (xc2.is_directed_right()) ?
        ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER) :
        ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER);
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 function object */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(this); }

  /*! A functor that compares the y-coordinates of two x-monotone arcs
   * immediately to the right of their intersection directional point.
   */
  class Compare_y_at_x_right_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_at_x_right_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Compare the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \param p the intersection point.
     * \return the relative position of xc1 with respect to xc2 immdiately to
     *         the right of p: SMALLER, EQUAL, or LARGER.
     * \pre the point p lies on both curves, and both of them must also be
     *      defined to its right (lexicographically).
     * \pre the arcs are not degenerate
     */
    Comparison_result operator()(const X_monotone_curve_2& xc1,
                                 const X_monotone_curve_2& xc2,
                                 const Point_2& p) const
    {
      CGAL_precondition(!xc1.is_degenerate());
      CGAL_precondition(!xc2.is_degenerate());

      // CGAL_precondition(p == xc1.left());
      // CGAL_precondition(p == xc2.left());

      // If Both arcs are vertical, they overlap:
      if (xc1.is_vertical() && xc2.is_vertical()) return EQUAL;
      if (xc1.is_vertical()) return LARGER;
      if (xc2.is_vertical()) return SMALLER;

      // Non of the arcs is verticel. Thus, non of the endpoints coincide with
      // a pole.
      // Compare the y-coord. at the x-coord of the most left right-endpoint.
      const Point_2& r1 = xc1.right();
      const Point_2& r2 = xc2.right();
      if (!r1.is_no_boundary()) {
        // use r2 and xc1:
        Oriented_side os = m_traits->oriented_side(xc1.normal(), r2);
        return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
          (xc1.is_directed_right()) ?
          ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER) :
          ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER);
      }
      if (!r2.is_no_boundary()) {
        // use r1 and xc2:
        Oriented_side os = m_traits->oriented_side(xc2.normal(), r1);
        return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
          (xc2.is_directed_right()) ?
          ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER) :
          ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER);
      }
      Comparison_result res = m_traits->compare_xy(r1, r2);
      if (res == LARGER) {
        // use r2 and xc1:
        Oriented_side os = m_traits->oriented_side(xc1.normal(), r2);
        return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
          (xc1.is_directed_right()) ?
          ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER) :
          ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER);
      }
      if (res == SMALLER) {
        // use r1 and xc2:
        Oriented_side os = m_traits->oriented_side(xc2.normal(), r1);
        return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
          (xc2.is_directed_right()) ?
          ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER) :
          ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER);
      }
      // res == equal
      // if p and r1 are antipodal, compare the plane normals
      const Kernel* kernel = m_traits;
      typename Kernel::Construct_opposite_direction_3 opposite_3 =
        kernel->construct_opposite_direction_3_object();
      Point_2 tmp1 = opposite_3(p);     // pacify msvc 10
      if (!kernel->equal_3_object()(Direction_3(tmp1), Direction_3(r1)))
        return EQUAL;

      Sign xsign = Traits::x_sign(p);
      Sign ysign = Traits::y_sign(p);
      Project project = (xsign == ZERO) ?
        ((ysign == POSITIVE) ? Traits::project_minus_xz : Traits::project_xz) :
        ((xsign == POSITIVE) ? Traits::project_yz : Traits::project_minus_yz);

      Direction_2 n1 = project(xc1.normal());
      Direction_2 n2 = project(xc2.normal());
      typename Kernel::Construct_opposite_direction_2 opposite_2 =
        kernel->construct_opposite_direction_2_object();
      if (!xc1.is_directed_right()) n1 = opposite_2(n1);
      if (!xc2.is_directed_right()) n2 = opposite_2(n2);
      if (kernel->equal_2_object()(n1, n2)) return EQUAL;
      const Direction_2 d(1, 0);
      return (kernel->counterclockwise_in_between_2_object()(n1, d, n2)) ?
        SMALLER : LARGER;
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 function object */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(this); }

  /*! A functor that checks whether two directional points and two x-monotone
   * arcs are identical.
   */
  class Equal_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Equal_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Determines whether the two x-monotone curves are the same (have the
     * same graph).
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \return true if the two curves are the same; false otherwise.
     */
    bool operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2) const
    {
      const Kernel* kernel = m_traits;
      typename Kernel::Equal_3 equal_3 = kernel->equal_3_object();
      if (xc1.is_full() || xc2.is_full()) {
        if (!xc1.is_full() || !xc2.is_full()) return false;
        typename Kernel::Construct_opposite_direction_3 opposite_3 =
          kernel->construct_opposite_direction_3_object();
        return (equal_3(xc1.normal(), xc2.normal()) ||
                equal_3(opposite_3(xc1.normal()), xc2.normal()));
      }

      if (xc1.is_meridian() || xc2.is_meridian()) {
        return (xc1.is_meridian() && xc2.is_meridian() &&
                equal_3(xc1.normal(), xc2.normal()));
      }

      return (equal_3(Direction_3(xc1.left()),  Direction_3(xc2.left())) &&
              equal_3(Direction_3(xc1.right()), Direction_3(xc2.right())));
    }

    /*! Determines whether the two points are the same.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return true if the two point are the same; false otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    {
      const Kernel* kernel = m_traits;
      return kernel->equal_3_object()(Direction_3(p1), Direction_3(p2));
    }
  };

  /*! Obtain an Equal_2 function object */
  Equal_2 equal_2_object() const { return Equal_2(this); }
  //@}

  /// \name Functor definitions to handle boundaries
  //@{

  /*! A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2 {
  public:
    /*! Obtains the parameter space at the end of an arc along the x-axis .
     * Note that if the arc end coincides with a pole, then unless the arc
     * coincides with the identification arc, the arc end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the arc coincides with the identification arc, it is assumed to
     * be smaller than any other object.
     * \param xcv the arc
     * \param ce the arc end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the arc xcv.
     *   ARR_LEFT_BOUNDARY  - the arc approaches the identification arc from
     *                        the right at the arc left end.
     *   ARR_INTERIOR       - the arc does not approache the identification arc.
     *   ARR_RIGHT_BOUNDARY - the arc approaches the identification arc from
     *                        the left at the arc right end.
     * \pre xcv does not coincide with the vertical identification arc.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const
    {
      if (xcv.is_vertical()) {
        CGAL_precondition(!xcv.is_on_boundary());
        return ARR_INTERIOR;
      }

      return (ce == ARR_MIN_END) ?
        ((xcv.left().is_no_boundary()) ? ARR_INTERIOR : ARR_LEFT_BOUNDARY) :
        ((xcv.right().is_no_boundary()) ? ARR_INTERIOR : ARR_RIGHT_BOUNDARY);
    }

    /*! Obtains the parameter space at a point along the x-axis.
     * \param p the point.
     * \return the parameter space at p.
     * \pre p does not lie on the vertical identification curve.
     */
    Arr_parameter_space operator()(const Point_2 p) const
    {
      /*! \todo For now the precondition is not applied. Instead, and as
       * a convention, if the point lies on the vertical identification
       * curve it is assumed to be smaller than anything else. Thus,
       * ARR_LEFT_BOUNDARY is returned.
       */
      return (p.is_mid_boundary()) ? ARR_LEFT_BOUNDARY : ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(); }

  /*! A function object that obtains the parameter space of a geometric
   * entity along the y-axis
   */
  class Parameter_space_in_y_2 {
  public:
    /*! Obtains the parameter space at the end of an arc along the y-axis .
     * Note that if the arc end coincides with a pole, then unless the arc
     * coincides with the identification arc, the arc end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the arc coincides with the identification arc, it is assumed to
     * be smaller than any other object.
     * \param xcv the arc
     * \param ce the arc end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the arc xcv.
     *   ARR_BOTTOM_BOUNDARY  - the arc approaches the south pole at the arc
     *                          left end.
     *   ARR_INTERIOR         - the arc does not approache a contraction point.
     *   ARR_TOP_BOUNDARY     - the arc approaches the north pole at the arc
     *                          right end.
     * There are no horizontal identification arcs!
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const
    {
      return (ce == ARR_MIN_END) ?
        ((xcv.left().is_min_boundary()) ?  ARR_BOTTOM_BOUNDARY: ARR_INTERIOR) :
        ((xcv.right().is_max_boundary()) ? ARR_TOP_BOUNDARY : ARR_INTERIOR);
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     * \pre p does not lie on the horizontal identification curve.
     * There are no horizontal identification arcs!
     */
    Arr_parameter_space operator()(const Point_2 p) const
    {
      return
        (p.is_min_boundary()) ? ARR_BOTTOM_BOUNDARY :
        (p.is_max_boundary()) ? ARR_TOP_BOUNDARY : ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }


  /*! A functor that compares the x-coordinate of arc ends and points on the
   * boundary of the parameter space.
   */
  class Compare_x_on_boundary_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_x_on_boundary_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Compare the x-limit of a direction with the x-coordinate of an
     * arc end on the boundary.
     * \param point the point direction.
     * \param xcv the arc, the endpoint of which is compared.
     * \param ce the arc-end indicator -
     *            ARR_MIN_END - the minimal end of xc or
     *            ARR_MAX_END - the maximal end of xc.
     * \return the comparison result:
     *         SMALLER - x(p) < x(xc, ce);
     *         EQUAL   - x(p) = x(xc, ce);
     *         LARGER  - x(p) > x(xc, ce).
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the arc xcv lies on a pole (implying ce is vertical).
     * \pre xcv does not coincide with the vertical identification curve.
     */
    Comparison_result operator()(const Point_2& point,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end CGAL_precondition_code(ce)) const
    {
      CGAL_precondition(point.is_no_boundary());
      CGAL_precondition_code
        (const Point_2& p2 = (ce == ARR_MIN_END) ? xcv.left() : xcv.right(););
      CGAL_precondition(!p2.is_no_boundary());
      CGAL_precondition(xcv.is_vertical());

      CGAL_precondition(!xcv.is_on_boundary());

      // xcv is vertical, but does not coincide with the discontinuity arc.
      // Obtain the direction contained in the underlying plane, which is
      // also on the xy-plane:
      Direction_3 normal = xcv.normal();
      Direction_2 q = (xcv.is_directed_right()) ?
        Direction_2(-(normal.dy()), normal.dx()) :
        Direction_2(normal.dy(), -(normal.dx()));
      Direction_2 p = Traits::project_xy(point);
      return m_traits->compare_x(p, q);
    }

    /*! Compare the x-coordinates of 2 arc ends near the boundary of the
     * parameter space.
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
     * \pre the ce1 end of the arc xcv1 lies on a pole (implying ce1 is
     *      vertical).
     * \pre the ce2 end of the arc xcv2 lies on a pole (implying ce2 is
     *      vertical).
     * \pre xcv1 does not coincide with the vertical identification curve.
     * \pre xcv2 does not coincide with the vertical identification curve.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end CGAL_precondition_code(ce1),
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end CGAL_precondition_code(ce2)) const
    {
      CGAL_precondition_code
        (const Point_2& p1 = (ce1 == ARR_MIN_END) ? xcv1.left() : xcv1.right(););
      CGAL_precondition(!p1.is_no_boundary());
      CGAL_precondition_code
        (const Point_2& p2 = (ce2 == ARR_MIN_END) ? xcv2.left() : xcv2.right(););
      CGAL_precondition(!p2.is_no_boundary());

      CGAL_precondition(xcv1.is_vertical());
      CGAL_precondition(xcv2.is_vertical());

      CGAL_precondition(!xcv1.is_on_boundary());
      CGAL_precondition(!xcv2.is_on_boundary());

      // Non of the arcs coincide with the discontinuity arc:
      // Obtain the directions contained in the underlying planes, which are
      // also on the xy-plane:
      Direction_3 normal1 = xcv1.normal();
      Direction_2 p = (xcv1.is_directed_right()) ?
        Direction_2(-(normal1.dy()), normal1.dx()) :
        Direction_2(normal1.dy(), -(normal1.dx()));
      Direction_3 normal2 = xcv2.normal();
      Direction_2 q = (xcv2.is_directed_right()) ?
        Direction_2(-(normal2.dy()), normal2.dx()) :
        Direction_2(normal2.dy(), -(normal2.dx()));
      return m_traits->compare_x(p, q);
    }

    /*! Compare the x-coordinate of two given points that lie on the
     * horizontal identification arc.
     * \param p1 the first point.
     * \param p2 the second point.
     * There is no horizontal identification arc!
     */
    Comparison_result operator()(const Point_2&, const Point_2&) const
    {
      CGAL_error_msg("There is no horizontal identification arc!");
      return SMALLER;
    }
  };

  /*! Obtain a Compare_x_on_boundary_2 function object */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
  { return Compare_x_on_boundary_2(this); }


  /*! A functor that compares the x-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_boundary_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_x_near_boundary_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:

    /*! Compare the x-coordinates of 2 arc ends near the boundary of the
     * parameter space.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the curve end indicator -
     *            ARR_MIN_END - the minimal end of curves or
     *            ARR_MAX_END - the maximal end of curves.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce) < x(xcv2, ce);
     *         EQUAL   - x(xcv1, ce) = x(xcv2, ce);
     *         LARGER  - x(xcv1, ce) > x(xcv2, ce).
     * \pre the ce end of the arc xcv1 lies on a pole.
     * \pre the ce end of the arc xcv2 lies on a pole.
     * \pre the the $x$-coordinates of xcv1 and xcv2 at their ce end are
     *      equal (implying that the curves overlap).
     * \pre xcv1 does not coincide with the vertical identification curve.
     * \pre xcv2 does not coincide with the vertical identification curve.
     */
    Comparison_result operator()(const X_monotone_curve_2&
                                   CGAL_precondition_code(xcv1),
                                 const X_monotone_curve_2&
                                   CGAL_precondition_code(xcv2),
                                 Arr_curve_end CGAL_precondition_code(ce)) const
    {
      CGAL_precondition_code
        (const Point_2& p1 = (ce == ARR_MIN_END) ? xcv1.left() : xcv1.right(););
      CGAL_precondition(!p1.is_no_boundary());
      CGAL_precondition_code
        (const Point_2& p2 = (ce == ARR_MIN_END) ? xcv2.left() : xcv2.right(););
      CGAL_precondition(!p2.is_no_boundary());

      CGAL_precondition(xcv1.is_vertical());
      CGAL_precondition(xcv2.is_vertical());
      CGAL_precondition(!xcv1.is_on_boundary());
      CGAL_precondition(!xcv2.is_on_boundary());

      return EQUAL;
    }
  };

  /*! Obtain a Compare_x_near_boundary_2 function object */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  { return Compare_x_near_boundary_2(this); }


  /*! A functor that compares the y-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_near_boundary_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Compare the y-coordinates of 2 curves at their ends near the boundary
     * of the parameter space.
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
      CGAL_precondition(!xcv1.is_degenerate());
      CGAL_precondition(!xcv2.is_degenerate());

      const Point_2& l1 = xcv1.left();
      const Point_2& r1 = xcv1.right();
      const Point_2& l2 = xcv2.left();
      const Point_2& r2 = xcv2.right();

      // If xcv1 is vertical, xcv1 coincides with the discontinuity arc:
      if (xcv1.is_vertical()) {
        CGAL_precondition(!l1.is_no_boundary());
        CGAL_precondition(!r1.is_no_boundary());
      }

      // If xcv2 is vertical, xcv2 coincides with the discontinuity arc:
      if (xcv2.is_vertical()) {
        CGAL_precondition(!l2.is_no_boundary());
        CGAL_precondition(!r2.is_no_boundary());
      }

      if (ce == ARR_MIN_END) {
        // Handle the south pole. It has the smallest y coords:
        if (l1.is_min_boundary())
          return (l2.is_min_boundary()) ? EQUAL : SMALLER;
        if (l2.is_min_boundary()) return LARGER;

        // None of xcv1 and xcv2 endpoints coincide with a pole:
        Comparison_result cr = m_traits->compare_y(l1, l2);
        if (cr != EQUAL) return cr;

        // If Both arcs are vertical, they overlap:
        if (xcv1.is_vertical() && xcv2.is_vertical()) return EQUAL;
        if (xcv1.is_vertical()) return LARGER;
        if (xcv2.is_vertical()) return SMALLER;

        // Non of the arcs is verticel. Thus, non of the endpoints coincide
        // with a pole.
        // Compare the y-coord. at the x-coord of the most left right-endpoint.
        CGAL_assertion(r1.is_no_boundary());
        CGAL_assertion(r2.is_no_boundary());

        if (m_traits->compare_xy(r1, r2) == LARGER) {
          // use r2 and xcv1:
          Oriented_side os =
            m_traits->oriented_side(xcv1.normal(), r2);
          return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
            (xcv1.is_directed_right()) ?
            ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER) :
            ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER);
        }
        // use r1 and xcv2:
        Oriented_side os = m_traits->oriented_side(xcv2.normal(), r1);
        return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
          (xcv2.is_directed_right()) ?
          ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER) :
          ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER);
      }

      // ce == ARR_MAX_END

      // Handle the north pole. It has the largest y coords:
      if (r1.is_max_boundary()) return (r2.is_max_boundary()) ? EQUAL : LARGER;
      if (r2.is_max_boundary()) return SMALLER;

      // None of xcv1 and xcv2 endpoints coincide with a pole:
      Direction_2 r1_xy = Traits::project_xy(r1);
      Comparison_result cr = m_traits->compare_y(r1, r2);
      if (cr != EQUAL) return cr;

      // If Both arcs are vertical, they overlap:
      if (xcv1.is_vertical() && xcv2.is_vertical()) return EQUAL;
      if (xcv1.is_vertical()) return LARGER;
      if (xcv2.is_vertical()) return SMALLER;

      // Compare to the left:
      Direction_2 p_r1 = Traits::project_xy(r1);
      cr = m_traits->compare_y(r1, r2);
      if (cr != EQUAL) return cr;

      // Non of the arcs is verticel. Thus, non of the endpoints coincide with
      // a pole.
      // Compare the y-coord. at the x-coord of the most right left-endpoint.
      CGAL_assertion(l1.is_no_boundary());
      CGAL_assertion(l2.is_no_boundary());

      if (m_traits->compare_xy(l1, l2) == SMALLER) {
        // use l2 and xcv1:
        Oriented_side os = m_traits->oriented_side(xcv1.normal(), l2);
        return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
          (xcv1.is_directed_right()) ?
          ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER) :
          ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER);
      }
      // use l1 and xcv2:
      Oriented_side os = m_traits->oriented_side(xcv2.normal(), l1);
      return (os == ON_ORIENTED_BOUNDARY) ? EQUAL :
        (xcv2.is_directed_right()) ?
        ((os == ON_NEGATIVE_SIDE) ? SMALLER : LARGER) :
        ((os == ON_NEGATIVE_SIDE) ? LARGER : SMALLER);
    }
  };

  /*! Obtain a Compare_y_near_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(this); }

  /*! A functor that indicates whether a geometric object lies on the
   * horizontal identification arc. In this setup there is no such entity.
   */
  class Is_on_x_identification_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

  public:
    /*! Determine whether a point lies on the horizontal identification arc.
     * \param p the point.
     * \return a Boolean indicating whether p lies on the horizontal
     * identification arc.
     */
    bool operator()(const Point_2&) const
    {
      CGAL_error_msg("There is no horizontal identification arc!");
      return false;
    }

    /*! Determine whether an arc coincides with the horizontal identification
     * arc.
     * \param xcv the arc.
     * \return a Boolean indicating whether xcv coincides with the horizontal
     * identification arc.
     */
    bool operator()(const X_monotone_curve_2&) const
    {
      CGAL_error_msg("There is no horizontal identification arc!");
      return false;
    }
  };

  /*! Obtain a Is_on_x_identification_2 function object */
  Is_on_x_identification_2 is_on_x_identification_2_object() const
  { return Is_on_x_identification_2(); }

  /*! A functor that indicates whether a geometric object lies on the
   * vertical identification arc.
   */
  class Is_on_y_identification_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Is_on_y_identification_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Determine whether a point lies on the vertical identification arc.
     * \param p the point.
     * \return a Boolean indicating whether p lies on the vertical
     * identification arc.
     */
    bool operator()(const Point_2& p) const { return p.is_mid_boundary(); }

    /*! Determine whether an arc coincides with the vertical identification
     * arc.
     * \param xcv the arc.
     * \return a Boolean indicating whether xcv coincides with the vertical
     * identification arc.
     */
    bool operator()(const X_monotone_curve_2& xcv) const
    { return xcv.is_on_boundary(); }
  };

  /*! Obtain a Is_on_y_identification_2 function object */
  Is_on_y_identification_2 is_on_y_identification_2_object() const
  { return Is_on_y_identification_2(this); }

  /*! A functor that compares the y-coordinate of two given points
   * that lie on the vertical identification arc.
   */
  class Compare_y_on_boundary_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Compare_y_on_boundary_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Compare the y-coordinate of two given points that lie on the vertical
     * identification curve.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return SMALLER - p1 is lexicographically smaller than p2;
     *         EQUAL   - p1 and p2 coincides;
     *         LARGER  - p1 is lexicographically larger than p2;
     * \pre p1 lies on the vertical identification arc.
     * \pre p2 lies on the vertical identification arc.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      CGAL_precondition(!p1.is_no_boundary());
      CGAL_precondition(!p2.is_no_boundary());
      return m_traits->compare_y(p1, p2);
    }
  };

  /*! Obtain a Compare_y_on_boundary_2 function object */
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const
  { return Compare_y_on_boundary_2(this); }
  //@}

  /// \name Functor definitions for supporting intersections.
  //@{

  /*! \class Make_x_monotone_2
   * A functor for subdividing arcs into x-monotone arcs that do not cross the
   * identification arc.
   */
  class Make_x_monotone_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Make_x_monotone_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Subdivide a given curve into x-monotone subcurves and insert them into
     * a given output iterator. As spherical_arcs are always x_monotone, only
     * one object will be contained in the iterator.
     * \param xc the curve.
     * \param oi the output iterator for the result. Its dereference type is a
     *           variant that wraps a \c Point_2 or an \c X_monotone_curve_2
     *           objects.
     * \return the past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& c, OutputIterator oi) const
    {
      typedef boost::variant<Point_2, X_monotone_curve_2>
        Make_x_monotone_result;

      if (c.is_degenerate()) {
        // The spherical_arc is a degenerate point - wrap it with an object:
        *oi++ = Make_x_monotone_result(c.right());
        return oi;
      }

      if (c.is_x_monotone()) {
        // The spherical arc is monotone - wrap it with an object:
        // *oi++ = make_object(X_monotone_curve_2(c));
        const X_monotone_curve_2* xc = &c;
        *oi++ = Make_x_monotone_result(*xc);
        return oi;
      }

      if (c.is_full()) {
        // The spherical arc is full
        if (c.is_vertical()) {
          // The arc is vertical => divide it into 2 meridians;
          const Direction_3& np = m_traits->neg_pole();
          const Direction_3& pp = m_traits->pos_pole();
          X_monotone_curve_2 xc1(np, pp, c.normal(), true, true);
          X_monotone_curve_2 xc2(pp, np, c.normal(), true, false);
          *oi++ = Make_x_monotone_result(xc1);
          *oi++ = Make_x_monotone_result(xc2);
          return oi;
        }
#if defined(CGAL_FULL_X_MONOTONE_GEODESIC_ARC_ON_SPHERE_IS_SUPPORTED)
        // The arc is not vertical => break it at the discontinuity arc:
        const X_monotone_curve_2 xc(c.normal());
        *oi++ = Make_x_monotone_result(xc);
#else
        // Full x-monotone arcs are not supported!
        // Split the arc at the intersection point with the complement of the
        // discontinuity arc:
        Direction_3 normal = c.normal();
        bool directed_right = Traits::x_sign(normal) == POSITIVE;
        Direction_3 d1(-(normal.dz()), 0, normal.dx());
        Direction_3 d2(normal.dz(), 0, -(normal.dx()));
        X_monotone_curve_2 xc1(d1, d2, normal, false, directed_right);
        X_monotone_curve_2 xc2(d2, d1, normal, false, directed_right);
        *oi++ = Make_x_monotone_result(xc1);
        *oi++ = Make_x_monotone_result(xc2);
#endif
        return oi;
      }

      const Point_2& source = c.source();
      const Point_2& target = c.target();
      const Direction_3& normal = c.normal();

      if (c.is_vertical()) {
        /* If one of the endpoints coincide with a pole, divide the arc at
         * the opposite pole:
         */
        const Direction_3& np = m_traits->neg_pole();
        const Direction_3& pp = m_traits->pos_pole();
        if (source.is_min_boundary() || target.is_min_boundary()) {
          X_monotone_curve_2 xc1(source, pp, normal, true, true);
          X_monotone_curve_2 xc2(pp, target, normal, true, false);
          *oi++ = Make_x_monotone_result(xc1);
          *oi++ = Make_x_monotone_result(xc2);
          return oi;
        }

        if (source.is_max_boundary() || target.is_max_boundary()) {
          X_monotone_curve_2 xc1(source, np, normal, true, false);
          X_monotone_curve_2 xc2(np, target, normal, true, true);
          *oi++ = Make_x_monotone_result(xc1);
          *oi++ = Make_x_monotone_result(xc2);
          return oi;
        }

        // None of the enpoints coincide with a pole.
        bool s_is_positive, t_is_positive, plane_is_positive;
        CGAL::Sign xsign = Traits::x_sign(normal);
        if (xsign == ZERO) {
          s_is_positive = Traits::x_sign(source) == POSITIVE;
          t_is_positive = Traits::x_sign(target) == POSITIVE;
          plane_is_positive = Traits::y_sign(normal) == NEGATIVE;
        } else {
          s_is_positive = Traits::y_sign(source) == POSITIVE;
          t_is_positive = Traits::y_sign(target) == POSITIVE;
          plane_is_positive = xsign == POSITIVE;
        }
        bool ccw = ((plane_is_positive && s_is_positive) ||
                    (!plane_is_positive && !s_is_positive));
        const Point_2& pole1 = (ccw) ? pp : np;
        X_monotone_curve_2 xc1(source, pole1, normal, true, ccw);
        *oi++ = Make_x_monotone_result(xc1);
        if (s_is_positive != t_is_positive) {
          // Construct 1 more arc:
          X_monotone_curve_2 xc2(pole1, target, normal, true, !ccw);
          *oi++ = Make_x_monotone_result(xc2);
          return oi;
        }
        // Construct 2 more arcs:
        const Point_2& pole2 = (ccw) ? np : pp;
        X_monotone_curve_2 xc2(pole1, pole2, normal, true, !ccw);
        *oi++ = Make_x_monotone_result(xc2);
        X_monotone_curve_2 xc3(pole2, target, normal, true, ccw);
        *oi++ = Make_x_monotone_result(xc3);
        return oi;
      }

      // The curve is not vertical, (none of the enpoints coincide with a pole)
#if (CGAL_IDENTIFICATION_XY == CGAL_X_MINUS_1_Y_0)
      Direction_3 dp = (CGAL::sign(normal.dz()) == POSITIVE) ?
        Direction_3(-(normal.dz()), 0, normal.dx()) :
        Direction_3(normal.dz(), 0, -(normal.dx()));
#else
      const Direction_2& xy = Traits::identification_xy();
      FT x = xy.dx();
      FT y = xy.dy();
      FT z((xy.dx() * normal.dx() + xy.dy() * normal.dy()) /  -(normal.dz()));
      Direction_3 dp(x, y, z);
#endif
      Point_2 p(dp, Point_2::MID_BOUNDARY_LOC);
      Direction_2 s = Traits::project_xy(source);
      Direction_2 t = Traits::project_xy(target);
      const Direction_2& d = Traits::identification_xy();
      const Kernel* kernel = m_traits;
      bool directed_right =
        kernel->counterclockwise_in_between_2_object()(d, s, t);

      X_monotone_curve_2 xc1(source, p, normal, false, directed_right);
      X_monotone_curve_2 xc2(p, target, normal, false, directed_right);
      *oi++ = Make_x_monotone_result(xc1);
      *oi++ = Make_x_monotone_result(xc2);
      return oi;
    }
  };

  /*! Obtain a Make_x_monotone_2 function object */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(this); }

  /*! A functor that splits an x-monotone arc at a directional point. */
  class Split_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Split_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Split a given x-monotone curve at a given point into two sub-curves.
     * \param xc the curve to split
     * \param p the split point.
     * \param xc1 (output) the left resulting subcurve. p is its right
     * endpoint.
     * \param xc2 (output) the right resulting subcurve. p is its left
     * endpoint.
     * \pre p lies on xc but is not one of its endpoints.
     * \pre xc is not degenerate
     */
    void operator()(const X_monotone_curve_2& xc, const Point_2& p,
                    X_monotone_curve_2& xc1, X_monotone_curve_2& xc2) const
    {
      CGAL_precondition(!xc.is_degenerate());
      const Point_2& source = xc.source();
      const Point_2& target = xc.target();
      CGAL_precondition_code(const Kernel* kernel = m_traits);
      CGAL_precondition_code
        (typename Kernel::Equal_3 equal_3 = kernel->equal_3_object());
      CGAL_precondition(!equal_3(Direction_3(p), Direction_3(source)));
      CGAL_precondition(!equal_3(Direction_3(p), Direction_3(target)));

      xc1.set_normal(xc.normal());
      xc1.set_is_vertical(xc.is_vertical());
      xc1.set_is_degenerate(false);
      xc1.set_is_empty(false);

      xc2.set_normal(xc.normal());
      xc2.set_is_vertical(xc.is_vertical());
      xc2.set_is_empty(false);

      if (xc.is_directed_right()) {
        xc1.set_source(source);
        xc1.set_target(p);
        xc1.set_is_directed_right(true);
        xc2.set_source(p);
        xc2.set_target(target);
        xc2.set_is_directed_right(true);
      } else {
        xc1.set_source(p);
        xc1.set_target(target);
        xc1.set_is_directed_right(false);
        xc2.set_source(source);
        xc2.set_target(p);
        xc2.set_is_directed_right(false);
      }
    }
  };

  /*! Obtain a Split_2 function object */
  Split_2 split_2_object() const { return Split_2(this); }

  /*! The clockwise-in-between function object */
  class Clockwise_in_between_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Clockwise_in_between_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    bool operator()(const Direction_2& d,
                    const Direction_2& d1, const Direction_2& d2) const
    {
      const Kernel* kernel = m_traits;
      return kernel->counterclockwise_in_between_2_object()(d, d2, d1);
    }
  };

  /*! Obtain a Clockwise_in_between function object */
  Clockwise_in_between_2 clockwise_in_between_2_object() const
  { return Clockwise_in_between_2(this); }

  /*! A functor that computes intersections between x-monotone arcs. */
  class Intersect_2 {
  private:

    /*! Computes the intersection between two arcs contained in the same plane
     * \param l1_3
     * \param r1_3
     * \param l2_3
     * \param r2_3
     * \param normal      - the normal of the common plane
     * \param vertical   - is the plane vertical
     * \param start      - the start 2d vertex
     * \param in_between - the in_between operator
     * \param project    - the projection function
     * \param oi         - the output iterator
     */
    template <typename In_between, typename OutputIterator>
    OutputIterator compute_intersection(const Point_2& l1_3,
                                        const Point_2 r1_3,
                                        const Point_2& l2_3,
                                        const Point_2 r2_3,
                                        const Direction_3& normal,
                                        bool vertical,
                                        const Direction_2& start,
                                        const In_between& in_between,
                                        Project project,
                                        OutputIterator oi) const
    {
      typedef std::pair<Point_2, Multiplicity>          Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;

      const Kernel* kernel = m_traits;
      typename Kernel::Equal_2 equal = kernel->equal_2_object();

      Direction_2 l1 = project(l1_3);
      Direction_2 r1 = project(r1_3);
      Direction_2 l2 = project(l2_3);
      Direction_2 r2 = project(r2_3);

      if (equal(l1, l2)) {
        const Point_2& trg = (in_between(r1, l2, r2)) ? r1_3 : r2_3;
        X_monotone_curve_2 xc(l1_3, trg, normal, vertical, true);
        *oi++ = Intersection_result(xc);
        return oi;
      }

      bool l1_eq_start = equal(l1, start);
      bool l2_eq_start = equal(l2, start);

      if (l1_eq_start || (!l2_eq_start && in_between(l1, start, l2))) {
        // The following applies only to full circles:
        if (l1_eq_start && equal(r2, start))
          *oi++ = Intersection_result(Intersection_point(r2_3, 1));
        if (in_between(r1, l1, l2)) return oi;      // no intersection
        if (equal(r1, l2)) {
          *oi++ = Intersection_result(Intersection_point(r1_3, 1));
          return oi;
        }
        const Point_2& trg = (in_between(r1, l2, r2)) ? r1_3 : r2_3;
        X_monotone_curve_2 xc(l2_3, trg, normal, vertical, true);
        *oi++ = Intersection_result(xc);
        return oi;
      }
      CGAL_assertion(l2_eq_start || in_between(l2, start, l1));
      // The following applies only to full circles:
      if (l2_eq_start && equal(r1, start))
        *oi++ = Intersection_result(Intersection_point(r1_3, 1));
      if (in_between(r2, l2, l1)) return oi;      // no intersection
      if (equal(r2, l1)) {
        *oi++ = Intersection_result(Intersection_point(r2_3, 1));
        return oi;
      }
      const Point_2& trg = (in_between(r1, l2, r2)) ? r1_3 : r2_3;
      X_monotone_curve_2 xc(l1_3, trg, normal, vertical, true);
      *oi++ = Intersection_result(xc);
      return oi;
    }

    /*! Determine whether a direction pierces an arc.
     * \param point the direction.
     * \param xc the arc.
     * \return true iff point pierces xc.
     * \pre point lies in the underlying plane of xc.
     */
    bool is_in_between(const Point_2& point,
                       const X_monotone_curve_2& xc) const
    {
      const Kernel* kernel = m_traits;
      CGAL_precondition(m_traits->has_on(xc.normal(), point));

      const Point_2& left = xc.left();
      const Point_2& right = xc.right();

      // Handle the poles:
      if (point.is_max_boundary()) return (right.is_max_boundary());
      if (point.is_min_boundary()) return (left.is_min_boundary());

      if (xc.is_vertical()) {
        // Compare the x coordinates. If they are not equal, return false:
        Direction_3 normal = xc.normal();
        bool plane_is_positive, p_is_positive;
        CGAL::Sign xsign = Traits::x_sign(normal);
        if (xsign == ZERO) {
          plane_is_positive = Traits::y_sign(normal) == NEGATIVE;
          p_is_positive = Traits::x_sign(point) == POSITIVE;
        } else {
          plane_is_positive = xsign == POSITIVE;
          p_is_positive = Traits::y_sign(point) == POSITIVE;
        }

        bool xc_is_positive = ((plane_is_positive && xc.is_directed_right()) ||
                               (!plane_is_positive && !xc.is_directed_right()));

        if ((xc_is_positive && !p_is_positive) ||
            (!xc_is_positive && p_is_positive))
          return false;

        // Compare the y-coords:
        return (((left.is_min_boundary()) ||
                 (m_traits->compare_y(point, left) != SMALLER)) &&
                ((right.is_max_boundary()) ||
                 (m_traits->compare_y(point, right) != LARGER)));
      }

      // The arc is not vertical. Compare the projections onto the xy-plane:
      typename Kernel::Equal_2 equal_2 = kernel->equal_2_object();
      Direction_2 p = Traits::project_xy(point);
      Direction_2 r = Traits::project_xy(right);
      if (equal_2(p, r)) return true;
      Direction_2 l = Traits::project_xy(left);
      if (equal_2(p, l)) return true;
      return kernel->counterclockwise_in_between_2_object()(p, l, r);
    }

  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Intersect_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Find the intersections of the two given curves and insert them into the
     * given output iterator. As two spherical_arcs may itersect only once,
     * only a single intersection will be contained in the iterator.
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \param oi the output iterator.
     * \return the past-the-end iterator.
     * \pre xc1 and xc2 are not degenerate
     */
    template<typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xc1,
                              const X_monotone_curve_2& xc2,
                              OutputIterator oi) const
    {
      CGAL_precondition(!xc1.is_degenerate());
      CGAL_precondition(!xc2.is_degenerate());

      typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;
      typedef typename Kernel::Counterclockwise_in_between_2
        Counterclockwise_in_between_2;
      typedef typename Kernel::Equal_3                  Equal_3;

      typedef std::pair<Point_2, Multiplicity>          Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;

      const Kernel* kernel = m_traits;

      Equal_3 equal_3 = kernel->equal_3_object();
      const Direction_3& normal1 = xc1.normal();
      const Direction_3& normal2 = xc2.normal();

      Direction_3 opposite_normal1 =
        kernel->construct_opposite_direction_3_object()(normal1);

      if (equal_3(normal1, normal2) || equal_3(opposite_normal1, normal2)) {
        // The underlying planes are the same
        Counterclockwise_in_between_2 ccib =
          kernel->counterclockwise_in_between_2_object();
        typename Traits::Clockwise_in_between_2 cib =
          m_traits->clockwise_in_between_2_object();

        if (xc1.is_vertical()) {
          // Both arcs are vertical
          bool res = kernel->equal_3_object()(normal1, normal2);
          if ((!res && (xc1.is_directed_right() == xc2.is_directed_right())) ||
              (res && (xc1.is_directed_right() != xc2.is_directed_right())))
          {
            if (xc1.left().is_min_boundary() && xc2.left().is_min_boundary())
              *oi++ = Intersection_result(Intersection_point(xc1.left(), 1));
            if (xc1.right().is_max_boundary() && xc2.right().is_max_boundary())
              *oi++ = Intersection_result(Intersection_point(xc1.right(), 1));
            return oi;
          }

          /*! If the endpoints of one arc coinside with the 2 poles resp,
           * the other arc is completely overlapping.
           */
          if (xc1.left().is_min_boundary() && xc1.right().is_max_boundary()) {
            *oi++ = Intersection_result(xc2);
            return oi;
          }
          if (xc2.left().is_min_boundary() && xc2.right().is_max_boundary()) {
            *oi++ = Intersection_result(xc1);
            return oi;
          }
          /*! Find an endpoint that does not coincide with a pole, and project
           * it onto the xy plane. If the projection coincide with the negative
           * x, project onto the zx plane. Otherwise project onto the yz plane.
           */
          const Point_2& point =
            xc1.left().is_min_boundary() ? xc1.right() : xc1.left();

          CGAL::Sign xsign = Traits::x_sign(normal1);
          bool xz_plane = xsign == ZERO;
          Project project =
            (xz_plane) ? Traits::project_xz : Traits::project_yz;

          Direction_3 normal = (xc1.is_directed_right()) ?
            normal1 : opposite_normal1;

          bool p_x_is_positive = Traits::x_sign(point) == POSITIVE;
          bool p_y_is_positive = Traits::y_sign(point) == POSITIVE;

          if ((xz_plane && p_x_is_positive) || (!xz_plane && p_y_is_positive)) {
            // The endpoints reside in the positive x-halfspace:
            return compute_intersection(xc1.left(), xc1.right(),
                                        xc2.left(), xc2.right(),
                                        normal, true, Traits::neg_y_2(),
                                        ccib, project, oi);
          }
          // The endpoints reside in the negative x-halfspace:
          return compute_intersection(xc1.left(), xc1.right(),
                                      xc2.left(), xc2.right(),
                                      normal, true, Traits::neg_y_2(),
                                      cib, project, oi);
        }

        // The arcs are not vertical:
        bool plane_is_positive = (Traits::z_sign(normal1) == POSITIVE);
        Direction_3 normal =
          (plane_is_positive) ? normal1 : opposite_normal1;
        return compute_intersection(xc1.left(), xc1.right(),
                                    xc2.left(), xc2.right(),
                                    normal, false, Traits::identification_xy(),
                                    ccib, Traits::project_xy, oi);
      }

      Vector_3 v =
        kernel->
        construct_cross_product_vector_3_object()(xc1.normal().vector(),
                                                  xc2.normal().vector());

      // Determine which one of the two directions:
      Point_2 ed(v.direction());
      if (is_in_between(ed, xc1) && is_in_between(ed, xc2)) {
        *oi++ = Intersection_result(Intersection_point(ed, 1));
        return oi;
      }

      Vector_3 vo(kernel->construct_opposite_vector_3_object()(v));
      Point_2 edo(vo.direction());
      if (is_in_between(edo, xc1) && is_in_between(edo, xc2)) {
        *oi++ = Intersection_result(Intersection_point(edo, 1));
        return oi;
      }
      return oi;
    }
  };

  /*! Obtain an Intersect_2 function object */
  Intersect_2 intersect_2_object() const { return Intersect_2(this); }

  /*! A functor that tests whether two x-monotone arcs can be merged. */
  class Are_mergeable_2 {
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Are_mergeable_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \return true if the two arcs are mergeable; false otherwise.
     * Two arcs are mergeable if:
     * 1. they are supported by the same plane, and
     * 2. share a common endpoint that is not on the identification arc
     */
    bool operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2) const
    {
      if (xc1.is_empty() || xc2.is_empty()) return true;
      if ((xc1.is_full() || xc1.is_meridian()) &&
          (xc2.is_full() || xc2.is_meridian())) return false;

      const Kernel* kernel = m_traits;
      typename Kernel::Equal_3 equal = kernel->equal_3_object();

      // Down cast to pass to kernel member functions
      const Direction_3& xc1_left = xc1.left();
      const Direction_3& xc2_left = xc2.left();
      const Direction_3& xc1_right = xc1.right();
      const Direction_3& xc2_right = xc2.right();
      if (xc1.is_degenerate() && xc2.is_degenerate())
        return equal(xc1_left, xc2_left);
      if ((xc1.is_full() || xc1.is_meridian()) && xc2.is_degenerate())
        return xc1.has_on(xc2.left());
      if ((xc2.is_full() || xc2.is_meridian()) && xc1.is_degenerate())
        return xc2.has_on(xc1.left());

      const Direction_3& normal1 = xc1.normal();
      const Direction_3& normal2 = xc2.normal();
      Direction_3 opposite_normal1 =
        kernel->construct_opposite_direction_3_object()(normal1);
      if (!equal(normal1, normal2) && !equal(opposite_normal1, normal2))
        return false;

      bool eq1 = equal(xc1_right, xc2_left);
      bool eq2 = equal(xc1_left, xc2_right);

#if defined(CGAL_FULL_X_MONOTONE_GEODESIC_ARC_ON_SPHERE_IS_SUPPORTED)
      if (eq1 && eq2) return true;
#else
      if (eq1 && eq2) return false;
#endif

      if (eq1 && xc2.left().is_no_boundary()) return true;
      if (eq2 && xc1.left().is_no_boundary()) return true;
      return false;
    }
  };

  /*! Obtain an Are_mergeable_2 function object */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(this); }

  /*! A functor that merges two x-monotone arcs into one */
  class Merge_2 {
  protected:
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_geodesic_arc_on_sphere_traits_2<Kernel>;

  public:
    /*! Merge two given x-monotone curves into a single curve (spherical_arc).
     * \param xc1 the first curve.
     * \param xc2 the second curve.
     * \param xc Output: the merged curve.
     * \pre the two curves are mergeable.
     */
    void operator()(const X_monotone_curve_2& xc1,
                    const X_monotone_curve_2& xc2,
                    X_monotone_curve_2& xc) const
    {
      CGAL_precondition (m_traits->are_mergeable_2_object()(xc1, xc2) == true);

      if (xc1.is_degenerate() || xc1.is_empty()) {
        xc = xc2;
        return;
      }

      if (xc2.is_degenerate() || xc2.is_empty()) {
        xc = xc1;
        return;
      }

      const Kernel* kernel = m_traits;
      typename Kernel::Equal_3 equal = kernel->equal_3_object();

      // Down cast to pass to kernel member functions
      const Direction_3& xc1_right = xc1.right();
      const Direction_3& xc2_left = xc2.left();

      xc.set_is_degenerate(false);
      xc.set_is_empty(false);
      xc.set_is_vertical(xc1.is_vertical());

      bool eq1 = equal(xc1_right, xc2_left);

#if defined(CGAL_FULL_X_MONOTONE_GEODESIC_ARC_ON_SPHERE_IS_SUPPORTED)
      const Direction_3& xc1_left = xc1.left();
      const Direction_3& xc2_right = xc2.right();
      bool eq2 = equal(xc1_left, xc2_right);
      if (eq1 && eq2) {
        const Point_2& p =
          xc1.source().is_mid_boundary() ? xc1.source() : xc1.target();
        xc.set_source(p);
        xc.set_target(p);
        xc.set_normal(xc1.normal());
        xc.set_is_full(true);
      }
#else
      CGAL_assertion_code(const Direction_3& xc1_left = xc1.left();
                          const Direction_3& xc2_right = xc2.right());
#endif
      if (xc1.is_directed_right() || xc2.is_directed_right()) {
        xc.set_normal(xc1.is_directed_right() ? xc1.normal() : xc2.normal());
        xc.set_is_directed_right(true);

        if (eq1) {
          xc.set_source(xc1.left());
          xc.set_target(xc2.right());
        }
        else {
          CGAL_assertion(equal(xc1_left, xc2_right));
          xc.set_source(xc2.left());
          xc.set_target(xc1.right());
        }
      }
      else {
        xc.set_normal(xc1.normal());
        xc.set_is_directed_right(false);

        if (eq1) {
          xc.set_source(xc2.right());
          xc.set_target(xc1.left());
        }
        else {
          CGAL_assertion(equal(xc1_left, xc2_right));
          xc.set_source(xc1.right());
          xc.set_target(xc2.left());
        }
      }
    }
  };

  /*! Obtain a Merge_2 function object */
  Merge_2 merge_2_object() const { return Merge_2(this); }
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

  /*! Obtain an Approximate_2 function object */
  Approximate_2 approximate_2_object() const { return Approximate_2(); }

  class Construct_x_monotone_curve_2 {
  public:

    /*! Return an x-monotone curve connecting the two given endpoints.
     * \param p the first point.
     * \param q the second point.
     * \pre p and q must not be the same.
     * \return a spherical_arc connecting p and q.
     */
    X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q) const
    { return X_monotone_curve_2(p, q); }
  };

  /*! Obtain a Construct_x_monotone_curve_2 function object */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(); }
  //@}


  /// \name Functor definitions for the Boolean set-operation traits.
  //@{

  class Compare_endpoints_xy_2 {
  public:

    /*!
     * Compare the endpoints of an $x$-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param xc the curve.
     * \return SMALLER if the curve is directed right;
     *         LARGER if the curve is directed left.
     */
    Comparison_result operator()(const X_monotone_curve_2& xc)
    { return (xc.is_directed_right()) ? SMALLER : LARGER; }
  };

  /*! Obtain a Compare_endpoints_xy_2 function object */
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

  /*! Obtain a Construct_opposite_2 function object */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }
  //@}

#if 0
  /*! Inserter for the spherical_arc class used by the traits-class */
  template <typename OutputStream>
  friend OutputStream& operator<<(OutputStream& os, const Point_2& p)
  {
    CGAL::To_double<typename Kernel::FT> todouble;
    os << static_cast<float>(todouble(p.dx())) << ", "
       << static_cast<float>(todouble(p.dy())) << ", "
       << static_cast<float>(todouble(p.dz()));
    return os;
  }

  /*! Inserter for the spherical_arc class used by the traits-class */
  template <typename OutputStream>
  friend OutputStream& operator<<(OutputStream& os,
                                  const X_monotone_curve_2& xc)
  {
    os << "(" << xc.left() << "), (" << xc.right() << ")";
    return os;
  }

  /*! Extractor for the spherical_arc class used by the traits-class */
  template <typename InputStream>
  friend InputStream& operator>>(InputStream& is, X_monotone_curve_2& arc)
  {
    CGAL_error_msg("Not implemented yet!");
    return is;
  }
#endif
};

/*! Represent an extended 3D direction that is used in turn to represent a
 * spherical-arc endpoint. The extended data consists of two flags that
 * indicate whether the point is on the x and on a y boundaries,
 * respectively.
 */
template <typename Kernel>
class Arr_extended_direction_3 : public Kernel::Direction_3 {
public:
  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::Direction_3                  Direction_3;

  /*! Enumeration of discontinuity type */
  enum Location_type {
    NO_BOUNDARY_LOC,
    MIN_BOUNDARY_LOC,
    MID_BOUNDARY_LOC,
    MAX_BOUNDARY_LOC
  };

private:
  typedef typename Kernel::Direction_2                  Direction_2;

  /*! The point discontinuity type */
  Location_type m_location;

  inline Sign x_sign(Direction_3 d) const { return CGAL::sign(d.dx()); }

  inline Sign y_sign(Direction_3 d) const { return CGAL::sign(d.dy()); }

  inline Sign z_sign(Direction_3 d) const { return CGAL::sign(d.dz()); }

public:
  /*! Default constructor */
  Arr_extended_direction_3() :
    Direction_3(0, 0, 1),
    m_location(MAX_BOUNDARY_LOC)
  {}

  /*! Constructor */
  Arr_extended_direction_3(const Direction_3& dir, Location_type location) :
    Direction_3(dir),
    m_location(location)
  {}

  /*! Constructor
   * \param x the x coordinate
   * \param y the y coordinate
   * \param z the z coordinate
   */
  Arr_extended_direction_3(const FT& x, const FT& y, const FT& z) :
    Direction_3(x, y, z)
  { init(); }

  /*! Constructor from a direction
   * \param dir the direction
   */
  Arr_extended_direction_3(const Direction_3& dir) : Direction_3(dir)
  { init(); }

  /*! Initialize from a direction
   * \param dir the direction
   */
  void init()
  {
    const Direction_3& dir = *this;
#if (CGAL_IDENTIFICATION_XY == CGAL_X_MINUS_1_Y_0)
    if (y_sign(dir) != ZERO) {
      m_location = NO_BOUNDARY_LOC;
      return;
    }
    CGAL::Sign signx = x_sign(dir);
    m_location =
      (signx == POSITIVE) ? NO_BOUNDARY_LOC :
      ((signx == NEGATIVE) ? MID_BOUNDARY_LOC :
       ((z_sign(dir) == NEGATIVE) ? MIN_BOUNDARY_LOC : MAX_BOUNDARY_LOC));
#else
    if ((x_sign(dir) == ZERO) && (y_sign(dir) == ZERO)) {
      m_location =
        (z_sign(dir) == NEGATIVE) ? MIN_BOUNDARY_LOC : MAX_BOUNDARY_LOC;
      return;
    }

    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;
    Direction_2 dir_xy = Traits::project_xy(dir);
    Kernel kernel;
    typename Kernel::Equal_2 equal_2 = kernel.equal_2_object();
    const Direction_2& xy = Traits::identification_xy();
    m_location = equal_2(dir_xy, xy) ? MID_BOUNDARY_LOC : NO_BOUNDARY_LOC;
#endif
  }

  /*! Copy constructor */
  Arr_extended_direction_3(const Arr_extended_direction_3& ed) :
    Direction_3(static_cast<const Direction_3&>(ed))
  {
    m_location = ed.discontinuity_type();
  }

  /*! Assignment operator */
  Arr_extended_direction_3& operator=(const Arr_extended_direction_3& ed)
  {
    *(static_cast<Direction_3*>(this)) = static_cast<const Direction_3&>(ed);
    m_location = ed.discontinuity_type();
    return (*this);
  }

  /*! Obtain the discontinuity type of the point */
  Location_type discontinuity_type() const
  { return m_location; }

  bool is_no_boundary() const { return (m_location == NO_BOUNDARY_LOC); }

  bool is_min_boundary() const { return (m_location == MIN_BOUNDARY_LOC); }

  bool is_mid_boundary() const { return (m_location == MID_BOUNDARY_LOC); }

  bool is_max_boundary() const { return (m_location == MAX_BOUNDARY_LOC); }
};

/*! A Representation of an x-monotone great circular arc embedded on a sphere,
 * as used by the Arr_geodesic_arc_on_sphere_traits_2 traits-class
 * An x-monotone great circular arc cannot cross the closed hemi-circle arc of
 * discontinuity, defined as the longitude that lies in the zx-plane, and is
 * contained in the open halfspace (x > 0).
 * \todo At this point such an arc cannot have an angle of 180 degrees.
 * \todo It is always directed from its source to its target.
 */
template <typename T_Kernel>
class Arr_x_monotone_geodesic_arc_on_sphere_3 {
public:
  typedef T_Kernel                                    Kernel;
  typedef typename Kernel::Direction_3                Direction_3;
  typedef typename Kernel::Plane_3                    Plane_3;
  typedef typename Kernel::Vector_3                   Vector_3;
  typedef typename Kernel::Direction_2                Direction_2;

protected:
  // For some reason compilation under Windows fails without the qualifier
  typedef CGAL::Arr_extended_direction_3<Kernel>      Arr_extended_direction_3;

  /*! The source point of the arc */
  Arr_extended_direction_3 m_source;

  /*! The target point of the arc */
  Arr_extended_direction_3 m_target;

  /*! The direction of the plane that contains the arc */
  Direction_3 m_normal;

  /*! The arc is vertical */
  bool m_is_vertical;

  /*! Target (lexicographically) larger than source */
  bool m_is_directed_right;

  /*! The arc is a full circle */
  bool m_is_full;

  /* The arc is degenerate - it consists of a single point */
  bool m_is_degenerate;

  /*! The arc is empty */
  bool m_is_empty;

  inline Sign x_sign(Direction_3 d) const { return CGAL::sign(d.dx()); }

  inline Sign y_sign(Direction_3 d) const { return CGAL::sign(d.dy()); }

  inline Sign z_sign(Direction_3 d) const { return CGAL::sign(d.dz()); }

  /*! Constructs a plane that contains two directions.
   * \param d1 the first direction.
   * \param d2 the second direction.
   */
  inline Direction_3 construct_normal_3(const Direction_3& d1,
                                        const Direction_3& d2) const
  {
    Kernel kernel;
    Vector_3 v = kernel.construct_cross_product_vector_3_object()(d1.vector(),
                                                                  d2.vector());
    return v.direction();
  }

public:
  /*! Default constructor - constructs an empty arc */
  Arr_x_monotone_geodesic_arc_on_sphere_3() :
    m_is_vertical(false),
    m_is_directed_right(false),
    m_is_full(false),
    m_is_degenerate(false),
    m_is_empty(true)
  {}

  /*! Constructor
   * \param src the source point of the arc
   * \param trg the target point of the arc
   * \param plane the plane that contains the arc
   * \param is_vertical is the arc vertical ?
   * \param is_directed_right is the arc directed from left to right?
   * \param is_full is the arc a full circle?
   * \param is_degenerate is the arc degenerate (single point)?
   */
  Arr_x_monotone_geodesic_arc_on_sphere_3
  (const Arr_extended_direction_3& src,
   const Arr_extended_direction_3& trg,
   const Direction_3& normal,
   bool is_vertical, bool is_directed_right,
   bool is_full = false, bool is_degenerate = false, bool is_empty = false) :
    m_source(src),
    m_target(trg),
    m_normal(normal),
    m_is_vertical(is_vertical),
    m_is_directed_right(is_directed_right),
    m_is_full(is_full),
    m_is_degenerate(is_degenerate),
    m_is_empty(is_empty)
  {}

  /*! Copy constructor
   * \param other the other arc
   */
  Arr_x_monotone_geodesic_arc_on_sphere_3
  (const Arr_x_monotone_geodesic_arc_on_sphere_3& other)
  {
    m_source = other.m_source;
    m_target = other.m_target;
    m_normal = other.m_normal;
    m_is_vertical = other.m_is_vertical;
    m_is_directed_right = other.m_is_directed_right;
    m_is_full = other.m_is_full;
    m_is_degenerate = other.m_is_degenerate;
    m_is_empty = other.m_is_empty;
  }

  /*! Assignment operator */
  Arr_x_monotone_geodesic_arc_on_sphere_3& operator=
  (const Arr_x_monotone_geodesic_arc_on_sphere_3& other)
  {
    m_source = other.m_source;
    m_target = other.m_target;
    m_normal = other.m_normal;
    m_is_vertical = other.m_is_vertical;
    m_is_directed_right = other.m_is_directed_right;
    m_is_full = other.m_is_full;
    m_is_degenerate = other.m_is_degenerate;
    m_is_empty = other.m_is_empty;
    return (*this);
  }

  /*! Construct the minor arc from two endpoint directions. The minor arc
   *  is the one with the smaller angle among the two geodesic arcs with
   * the given endpoints.
   * 1. Find out whether the arc is x-monotone.
   * 2. If it is x-monotone,
   *    2.1 Find out whether it is vertical, and
   *    2.2 whether the target is larger than the source (directed right).
   * The arc is vertical, iff
   * 1. one of its endpoint direction pierces a pole, or
   * 2. the projections onto the xy-plane coincide.
   * \param source the source point.
   * \param target the target point.
   * \pre the source and target cannot be equal.
   * \pre the source and target cannot be antipodal.
   */
  Arr_x_monotone_geodesic_arc_on_sphere_3
  (const Arr_extended_direction_3& source,
   const Arr_extended_direction_3& target) :
    m_source(source),
    m_target(target),
    m_is_full(false),
    m_is_degenerate(false),
    m_is_empty(false)
  {
    // MSVC 10 complains when the casting below is not present probably due
    // to a bug (in MSVC 10).
    CGAL_precondition_code(Kernel kernel);
    CGAL_precondition(!kernel.equal_3_object()
                      (kernel.construct_opposite_direction_3_object()(m_source),
                       (const typename Kernel::Direction_3&)(m_target)));
    m_normal = construct_normal_3(m_source, m_target);
    init();
  }

  /*! Initialize a spherical_arc given that the two endpoint directions
   * have been set. It is assumed that the arc is the one with the smaller
   * angle among the two.
   * 1. Find out whether the arc is x-monotone.
   * 2. If it is x-monotone,
   *    2.1 Find out whether it is vertical, and
   *    2.2 whether the target is larger than the source (directed right).
   * The arc is vertical, iff
   * 1. one of its endpoint direction pierces a pole, or
   * 2. the projections onto the xy-plane coincide.
   * \param source the source point.
   * \param target the target point.
   * \pre the source and target cannot be equal.
   */
  void init()
  {
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    Kernel kernel;
    CGAL_precondition(!kernel.equal_3_object()(Direction_3(m_source),
                                               Direction_3(m_target)));

    // Check whether any one of the endpoint coincide with a pole:
    if (m_source.is_max_boundary()) {
      set_is_vertical(true);
      set_is_directed_right(false);
      return;
    }
    if (m_source.is_min_boundary()) {
      set_is_vertical(true);
      set_is_directed_right(true);
      return;
    }
    if (m_target.is_max_boundary()) {
      set_is_vertical(true);
      set_is_directed_right(true);
      return;
    }
    if (m_target.is_min_boundary()) {
      set_is_vertical(true);
      set_is_directed_right(false);
      return;
    }

    // None of the enpoints coincide with a pole:
    Direction_2 s = Traits::project_xy(m_source);
    Direction_2 t = Traits::project_xy(m_target);

    Orientation orient = Traits::orientation(s, t);
    if (orient == COLLINEAR) {
      set_is_vertical(true);
      const Direction_2& nx = Traits::neg_x_2();
      if (Traits::orientation(nx, s) == COLLINEAR) {
        // Project onto xz plane:
        s = Traits::project_xz(m_source);
        t = Traits::project_xz(m_target);
        const Direction_2& ny = Traits::neg_y_2();
        Orientation orient1 = Traits::orientation(ny, s);
        CGAL_assertion_code(Orientation orient2 = Traits::orientation(ny, t));
        CGAL_assertion(orient1 == orient2);
        orient = Traits::orientation(s, t);
        CGAL_assertion(orient != COLLINEAR);
        if (orient1 == LEFT_TURN) {
          set_is_directed_right(orient == LEFT_TURN);
          return;
        }
        set_is_directed_right(orient == RIGHT_TURN);
        return;
      }
      // Project onto yz plane:
      s = Traits::project_yz(m_source);
      t = Traits::project_yz(m_target);
      const Direction_2& ny = Traits::neg_y_2();
      Orientation orient1 = Traits::orientation(ny, s);
      CGAL_assertion_code(Orientation orient2 = Traits::orientation(ny, t));
      CGAL_assertion(orient1 == orient2);
      if (orient1 == LEFT_TURN) {
        orient = Traits::orientation(s, t);
        CGAL_assertion(orient != COLLINEAR);
        set_is_directed_right(orient == LEFT_TURN);
        return;
      }
      orient = Traits::orientation(s, t);
      CGAL_assertion(orient != COLLINEAR);
      set_is_directed_right(orient == RIGHT_TURN);
      return;
    }

    // The arc is not vertical!
    set_is_vertical(false);
    set_is_directed_right(orient == LEFT_TURN);
    set_is_full(kernel.equal_3_object()(Direction_3(m_source),
                                        Direction_3(m_target)));
  }

  /*! Construct a full spherical_arc from a plane
   * \param plane the containing plane.
   * \pre the plane is not vertical
   */
  Arr_x_monotone_geodesic_arc_on_sphere_3(const Direction_3& normal) :
    m_normal(normal),
    m_is_vertical(false),
    m_is_directed_right(z_sign(normal) == POSITIVE),
    m_is_full(true),
    m_is_degenerate(false),
    m_is_empty(false)
  {
    CGAL_precondition(z_sign(normal) != ZERO);

#if (CGAL_IDENTIFICATION_XY == CGAL_X_MINUS_1_Y_0)
    Direction_3 d = (CGAL::sign(normal.dz()) == POSITIVE) ?
      Direction_3(-(normal.dz()), 0, normal.dx()) :
      Direction_3(normal.dz(), 0, -(normal.dx()));
#else
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;
    typedef typename Kernel::FT                         FT;

    const Direction_2& xy = Traits::identification_xy();
    FT x = xy.dx();
    FT y = xy.dy();
    FT z((xy.dx() * normal.dx() + xy.dy() * normal.dy()) / -(normal.dz()));
    Direction_3 d(x, y, z);
#endif
    m_source = m_target =
      Arr_extended_direction_3(d, Arr_extended_direction_3::MID_BOUNDARY_LOC);
  }

  /*! Construct a full spherical_arc from a common endpoint and a plane
   * \param plane the containing plane.
   * \pre the point lies on the plane
   * \pre the point lies on the open discontinuity arc
   */
  Arr_x_monotone_geodesic_arc_on_sphere_3
  (const Arr_extended_direction_3& point,
   const Direction_3& normal) :
    m_source(point),
    m_target(point),
    m_normal(normal),
    m_is_vertical(false),
    m_is_directed_right(z_sign(normal) == POSITIVE),
    m_is_full(true),
    m_is_degenerate(false),
    m_is_empty(false)
  {
    CGAL_precondition(has_on(point));
    CGAL_precondition(z_sign(normal) != ZERO);
#if !defined(CGAL_FULL_X_MONOTONE_GEODESIC_ARC_ON_SPHERE_IS_SUPPORTED)
    CGAL_error_msg( "Full x-monotone arcs are not supported!");
#endif
  }

  /*! Construct a spherical_arc from two endpoints directions contained
   * in a plane.
   * \param plane the containing plane.
   * \param source the source-point direction.
   * \param target the target-point direction.
   * \pre Both endpoint lie on the given plane.
   * \pre Both endpoint lie on the given plane.
   */
  Arr_x_monotone_geodesic_arc_on_sphere_3
  (const Arr_extended_direction_3& source,
   const Arr_extended_direction_3& target,
   const Direction_3& normal) :
    m_source(source),
    m_target(target),
    m_normal(normal),
    m_is_full(false),
    m_is_degenerate(false),
    m_is_empty(false)
  {
    CGAL_precondition(has_on(source));
    CGAL_precondition(has_on(target));

    // Check whether any one of the endpoint coincide with a pole:
    if (source.is_max_boundary()) {
      set_is_vertical(true);
      set_is_directed_right(false);
      return;
    }
    if (source.is_min_boundary()) {
      set_is_vertical(true);
      set_is_directed_right(true);
      return;
    }
    if (target.is_max_boundary()) {
      set_is_vertical(true);
      set_is_directed_right(true);
      return;
    }
    if (target.is_min_boundary()) {
      set_is_vertical(true);
      set_is_directed_right(false);
      return;
    }

    if (z_sign(normal) == ZERO) {
      set_is_vertical(true);
      bool s_is_positive, plane_is_positive;
      CGAL::Sign xsign = x_sign(normal);
      if (xsign == ZERO) {
        s_is_positive = x_sign(source) == POSITIVE;
        plane_is_positive = y_sign(normal) == NEGATIVE;
      } else {
        s_is_positive = y_sign(source) == POSITIVE;
        plane_is_positive = xsign == POSITIVE;
      }
      bool ccw = ((plane_is_positive && s_is_positive) ||
                  (!plane_is_positive && !s_is_positive));
      set_is_directed_right(ccw);
      return;
    }

    // The arc is not vertical!
    set_is_vertical(false);
    set_is_directed_right(z_sign(normal) == POSITIVE);
  }

  /*! Set the source endpoint direction.
   * \param p the endpoint to set.
   */
  void set_source(const Direction_3& p) { m_source = p; }

  /*! Set the target endpoint direction.
   * \param p the endpoint to set.
   */
  void set_target(const Direction_3& p) { m_target = p; }

  /*! Set the direction of the underlying plane.
   * \param normal the plane direction.
   */
  void set_normal(const Direction_3& normal) { m_normal = normal; }

  void set_is_vertical(bool flag) { m_is_vertical = flag; }
  void set_is_directed_right(bool flag) { m_is_directed_right = flag; }
  void set_is_full(bool flag) { m_is_full = flag; }
  void set_is_degenerate(bool flag) { m_is_degenerate = flag; }
  void set_is_empty(bool flag) { m_is_empty = flag; }

  /*! Obtain the source */
  const Arr_extended_direction_3& source() const { return m_source; }

  /*! Obtain the target */
  const Arr_extended_direction_3& target() const { return m_target; }

  /*! Obtain the normal to the containing plane */
  const Direction_3& normal() const { return m_normal; }

  /*! Obtain the (lexicographically) left endpoint direction */
  const Arr_extended_direction_3& left() const
  { return (m_is_directed_right ? m_source : m_target); }

  /*! Obtain the (lexicographically) right endpoint */
  const Arr_extended_direction_3& right() const
  { return (m_is_directed_right ? m_target : m_source); }

  /*! Determines whether the curve is vertical */
  bool is_vertical() const { return m_is_vertical; }

  /*! Determines whether the curve is directed lexicographically from left to
   * right
   */
  bool is_directed_right() const { return m_is_directed_right; }

  /*! Determines whether the curve is a full circle */
  bool is_full() const { return m_is_full; }

  /*! Determines whether the curve is degenerate */
  bool is_degenerate() const { return m_is_degenerate; }

  /*! Determines whether the curve is degenerate */
  bool is_empty() const { return m_is_empty; }

  /*! Determine whether the curve lie on the identification arc */
  bool is_on_boundary() const
  {
    /* If the curve is not vertical and non of its endpoints lie on the
     * boundary, the arc itself cannot lie on the identification arc.
     */
    if (m_source.is_no_boundary() || m_target.is_no_boundary() ||
        !is_vertical())
      return false;

    /*! The curve is vertical. If at least one endpoint lies on the open
     * identification arc, it entirely lies on it.
     */
    if (m_source.is_mid_boundary() || m_target.is_mid_boundary())
      return true;

    /* Both endpoints lie on opposite poles respectively. If the normal
     * coincides with the normal of the plane that contains the identification
     * arc, the arc lies on the identification arc.
     */
#if (CGAL_IDENTIFICATION_XY == CGAL_X_MINUS_1_Y_0)
    return ((x_sign(m_normal) == ZERO) &&
            (((y_sign(m_normal) == NEGATIVE) && !is_directed_right()) ||
             ((y_sign(m_normal) == POSITIVE) && is_directed_right())));
#else
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    const Direction_3& iden_normal = Traits::identification_normal();
    const Direction_2 iden_normal_xy = Traits::project_xy(iden_normal);
    Direction_2 normal_xy = Traits::project_xy(m_normal);
    Kernel kernel;
    if (is_directed_right()) {
      return kernel.equal_2_object()(normal_xy, iden_normal_xy);
    } else {
      Direction_2 opposite_normal_xy =
        kernel.construct_opposite_direction_2_object()(normal_xy);
      return kernel.equal_2_object()(opposite_normal_xy, iden_normal_xy);
    }
#endif
  }

  /*! Determine whether the given point is in the x-range of the
   * spherical_arc.
   * \param point the query point direction.
   * \return true if point is in the x-range of the (closed) spherical_arc and
   * false otherwise.
   * \pre point does not coincide with one of the poles
   */
  bool is_in_x_range(const Arr_extended_direction_3& point) const
  {
    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    CGAL_precondition(!point.is_min_boundary());
    CGAL_precondition(!point.is_max_boundary());

    Direction_2 p = Traits::project_xy(point);
    if (is_vertical()) {
      Direction_2 q = (is_directed_right()) ?
        Direction_2(-(m_normal.dy()), m_normal.dx()) :
        Direction_2(m_normal.dy(), -(m_normal.dx()));
      Kernel kernel;
      return kernel.equal_2_object()(p, q);
    }

    // The curve is not vertical:
    Direction_2 r = Traits::project_xy(right());
    Kernel kernel;
    if (kernel.equal_2_object()(p, r)) return true;
    Direction_2 l = Traits::project_xy(left());
    if (kernel.equal_2_object()(p, l)) return true;
    return kernel.counterclockwise_in_between_2_object()(p, l, r);
  }

  /*! Determines whether the curve is a meridian */
  bool is_meridian() const
  { return left().is_min_boundary() && right().is_max_boundary(); }

#if 0
  /*! Create a bounding box for the spherical_arc */
  Bbox_2 bbox() const
  {
    Kernel kernel;
    Segment_2 seg = kernel.construct_spherical_arc_2_object()(this->m_source,
                                                              this->m_target);
    return kernel.construct_bbox_2_object()(seg);
  }
#endif

  /*! Flip the spherical_arc (swap it source and target) */
  Arr_x_monotone_geodesic_arc_on_sphere_3 opposite() const
  {
    Arr_x_monotone_geodesic_arc_on_sphere_3 opp;
    opp.m_source = this->m_target;
    opp.m_target = this->m_source;
    opp.m_normal = this->m_normal;
    opp.m_is_directed_right = !(this->is_directed_right());
    opp.m_is_vertical = this->is_vertical();
    opp.m_is_full = this->is_full();
    opp.m_is_degenerate = this->is_degenerate();
    opp.m_is_empty = this->is_empty();
    return opp;
  }

  /*! Determined whether a direction is contained in a plane
   * \param plane the 3D plane.
   * \param dir the 3D direction.
   * \return true if dir is contained in plane; false otherwise.
   * \pre the plane contains the origin.
   */
  inline bool has_on(const Direction_3& dir) const
  {
    typename Kernel::FT dot = normal().vector() * dir.vector();
    return CGAL::sign(dot) == ZERO;
  }
};

/*! A representation of a general great circular arc embedded on a sphere,
 * used by the Arr_geodesic_arc_on_sphere_traits_2 traits-class
 * An arc is uniqely represented by a plane p, and two endpoints the source
 * s and the target t, which lie in the plane p. The points of the arc are
 * the locus of points visited when moving from the source s toward the
 * target t on the plane p in counterclockwise direction along the circle
 * defined by s and t.
 */
template <typename T_Kernel>
class Arr_geodesic_arc_on_sphere_3 :
  public Arr_x_monotone_geodesic_arc_on_sphere_3<T_Kernel> {
public:
  typedef T_Kernel                                              Kernel;

protected:
  typedef Arr_x_monotone_geodesic_arc_on_sphere_3<Kernel>       Base;

public:
  typedef typename Base::Plane_3                                Plane_3;
  typedef typename Base::Direction_3                            Direction_3;
  typedef typename Base::Direction_2                            Direction_2;

protected:
  // For some reason compilation under Windows fails without the qualifier
  typedef CGAL::Arr_extended_direction_3<Kernel>    Arr_extended_direction_3;

  using Base::x_sign;
  using Base::y_sign;
  using Base::z_sign;

  /*! Indicates whether the arc is x-monotone */
  bool m_is_x_monotone;

public:
  /*! Default constructor - constructs an empty arc */
  Arr_geodesic_arc_on_sphere_3() : Base(), m_is_x_monotone(true) {}

  /*! Copy constructor
   * \param other the other arc
   */
#ifdef DOXYGEN_RUNNING
  Arr_geodesic_arc_on_sphere_3
  (const Arr_geodesic_arc_on_sphere_3& other) : Base(other)
  {
    m_is_x_monotone = other.m_is_x_monotone;
  }
#endif

  /*! Constructor
   * \param src the source point of the arc
   * \param trg the target point of the arc
   * \param normal the normal to the  plane that contains the arc
   * \param is_x_monotone is arc  x-monotone ?
   * \param is_vertical is the arc vertical ?
   * \param is_directed_right is the arc directed from left to right?
   * \param is_full is the arc a full (great) circle?
   * \param is_degenerate is the arc degenerate (single point)?
   * \pre plane contains the origin
   * \pre plane contains src
   * \pre plane contains trg
   */
  Arr_geodesic_arc_on_sphere_3(const Arr_extended_direction_3& src,
                               const Arr_extended_direction_3& trg,
                               const Direction_3& normal,
                               bool is_x_monotone, bool is_vertical,
                               bool is_directed_right,
                               bool is_full = false,
                               bool is_degenerate = false,
                               bool is_empty = false) :
    Base(src, trg, normal,
         is_vertical, is_directed_right, is_full, is_degenerate, is_empty),
    m_is_x_monotone(is_x_monotone)
  {
    CGAL_precondition(this->has_on(src));
    CGAL_precondition(this->has_on(trg));
  }

  /*! Construct a spherical_arc from two endpoint directions. It is assumed
   * that the arc is the one with the smaller angle among the two.
   * 1. Find out whether the arc is x-monotone.
   * 2. If it is x-monotone,
   *    2.1 Find out whether it is vertical, and
   *    2.2 whether the target is larger than the source (directed right).
   * The arc is vertical, iff
   * 1. one of its endpoint direction pierces a pole, or
   * 2. the projections onto the xy-plane coincide.
   * \param source the source point.
   * \param target the target point.
   * \pre the source and target cannot be equal.
   * \pre the source and target cannot be the opoosite of each other.
   */
  Arr_geodesic_arc_on_sphere_3(const Arr_extended_direction_3& source,
                               const Arr_extended_direction_3& target)
  {
    this->set_source(source);
    this->set_target(target);
    this->set_is_full(false);
    this->set_is_degenerate(false);
    this->set_is_empty(false);

    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel>         Traits;
    typedef typename Kernel::Direction_2                        Direction_2;
    typedef typename Kernel::Direction_3                        Direction_3;

    Kernel kernel;
    CGAL_precondition(!kernel.equal_3_object()(Direction_3(source),
                                               Direction_3(target)));
    CGAL_precondition(!kernel.equal_3_object()
                      (kernel.construct_opposite_direction_3_object()(source),
                       static_cast<const Direction_3&>(target)));
    this->m_normal = this->construct_normal_3(source, target);

    // Check whether one of the endpoints coincides with a pole: */
    if (source.is_max_boundary()) {
      this->set_is_vertical(true);
      this->set_is_directed_right(false);
      set_is_x_monotone(true);
      return;
    }
    if (source.is_min_boundary()) {
      this->set_is_vertical(true);
      this->set_is_directed_right(true);
      set_is_x_monotone(true);
      return;
    }
    if (target.is_max_boundary()) {
      this->set_is_vertical(true);
      this->set_is_directed_right(true);
      set_is_x_monotone(true);
      return;
    }
    if (target.is_min_boundary()) {
      this->set_is_vertical(true);
      this->set_is_directed_right(false);
      set_is_x_monotone(true);
      return;
    }

    // None of the enpoints coincide with a pole:
    Direction_3 normal = this->m_normal;
    if (z_sign(normal) == ZERO) {
      // The arc is vertical
      this->set_is_vertical(true);
      bool s_is_positive, t_is_positive, plane_is_positive;
      CGAL::Sign xsign = x_sign(normal);
      if (xsign == ZERO) {
        s_is_positive = x_sign(source) == POSITIVE;
        t_is_positive = x_sign(target) == POSITIVE;
        plane_is_positive = y_sign(normal) == NEGATIVE;
      } else {
        s_is_positive = y_sign(source) == POSITIVE;
        t_is_positive = y_sign(target) == POSITIVE;
        plane_is_positive = xsign == POSITIVE;
      }
      set_is_x_monotone(s_is_positive == t_is_positive);
      bool ccw = ((plane_is_positive && s_is_positive) ||
                  (!plane_is_positive && !s_is_positive));
      this->set_is_directed_right(ccw);
      return;
    }

    // The arc is not vertical!
    this->set_is_vertical(false);
    Direction_2 s = Traits::project_xy(source);
    Direction_2 t = Traits::project_xy(target);
    Orientation orient = Traits::orientation(s, t);

    const Direction_2& d = Traits::identification_xy();
    if (orient == LEFT_TURN) {
      this->set_is_directed_right(true);
      set_is_x_monotone(!kernel.counterclockwise_in_between_2_object()(d, s, t));
      return;
    }

    // (orient == RIGHT_TURN)
    this->set_is_directed_right(false);
    set_is_x_monotone(!kernel.counterclockwise_in_between_2_object()(d, t, s));
    return;
  }

  /*! Construct a spherical_arc from two endpoint directions contained
   * in a plane.
   * \param plane the containing plane.
   * \param source the source-point direction.
   * \param target the target-point direction.
   * \param normal the normal to the plane containing the arc
   * \pre plane contain the origin
   * \pre Both endpoints lie on the given plane.
   */
  Arr_geodesic_arc_on_sphere_3(const Arr_extended_direction_3& source,
                               const Arr_extended_direction_3& target,
                               const Direction_3& normal)
  {
    Kernel kernel;

    this->set_source(source);
    this->set_target(target);
    this->set_normal(normal);
    this->set_is_degenerate(false);
    this->set_is_empty(false);

    typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel> Traits;

    CGAL_precondition(this->has_on(source));
    CGAL_precondition(this->has_on(target));

    if (z_sign(normal) == ZERO) {
      this->set_is_vertical(true);

      // Check whether both endpoint coincide with the poles:
      if (source.is_min_boundary() && target.is_max_boundary()) {
        // Both endpoints coincide with the 2 poles respectively.
        this->set_is_directed_right(true);
        this->set_is_full(false);
        set_is_x_monotone(true);
        return;
      }

      if (source.is_max_boundary() && target.is_min_boundary()) {
        // Both endpoints coincide with the 2 poles respectively.
        this->set_is_directed_right(false);
        this->set_is_full(false);
        set_is_x_monotone(true);
        return;
      }

      CGAL::Sign xsign = x_sign(normal);
      bool xz_plane = xsign == ZERO;
      bool s_is_positive, t_is_positive, plane_is_positive;
      if (xz_plane) {
        s_is_positive = x_sign(source) == POSITIVE;
        t_is_positive = x_sign(target) == POSITIVE;
        plane_is_positive = y_sign(normal) == NEGATIVE;
      } else {
        s_is_positive = y_sign(source) == POSITIVE;
        t_is_positive = y_sign(target) == POSITIVE;
        plane_is_positive = xsign == POSITIVE;
      }

      // Process degenerate cases:
      if (source.is_min_boundary()) {
        this->set_is_directed_right(true);
        set_is_x_monotone((plane_is_positive && t_is_positive) ||
                          (!plane_is_positive && !t_is_positive));
        return;
      }
      if (source.is_max_boundary()) {
        this->set_is_directed_right(false);
        set_is_x_monotone((plane_is_positive && !t_is_positive) ||
                          (!plane_is_positive && t_is_positive));
        return;
      }
      if (target.is_min_boundary()) {
        this->set_is_directed_right(false);
        set_is_x_monotone((plane_is_positive && !s_is_positive) ||
                          (!plane_is_positive && s_is_positive));
        return;
      }
      if (target.is_max_boundary()) {
        this->set_is_directed_right(true);
        set_is_x_monotone((plane_is_positive && s_is_positive) ||
                          (!plane_is_positive && !s_is_positive));
        return;
      }
      if (s_is_positive != t_is_positive) {
        set_is_x_monotone(false);
        return;
      }

      /* Non of the endpoints coincide with a pole.
       * The projections of both endpoints lie on the same hemi-circle.
       * Thus, either the arc is x-monotone, or it includes both poles.
       * This means that it is sufficient to check whether one pole lies
       * on the arc in order to determine x-monotonicity
       */

      typename Traits::Project project =
        (xz_plane) ? Traits::project_xz : Traits::project_yz;
      Direction_2 s = project(source);
      Direction_2 t = project(target);
      const Direction_2& ny = Traits::neg_y_2();
      typename Kernel::Counterclockwise_in_between_2 ccib =
        kernel.counterclockwise_in_between_2_object();
      set_is_x_monotone((plane_is_positive && !ccib(ny, s, t)) ||
                        (!plane_is_positive && !ccib(ny, t, s)));

      bool ccw = ((plane_is_positive && s_is_positive) ||
                  (!plane_is_positive && !s_is_positive));
      this->set_is_directed_right(ccw);
      return;
    }

    // The arc is not vertical!
    this->set_is_vertical(false);
    this->set_is_directed_right(z_sign(normal) == POSITIVE);
    const Direction_2& d = Traits::identification_xy();
    Direction_2 s = Traits::project_xy(source);
    Direction_2 t = Traits::project_xy(target);
    typename Kernel::Counterclockwise_in_between_2 ccib =
      kernel.counterclockwise_in_between_2_object();
    bool plane_is_positive = (z_sign(normal) == POSITIVE);
    set_is_x_monotone((plane_is_positive && !ccib(d, s, t)) ||
                      (!plane_is_positive && !ccib(d, t, s)));
  }

  /*! Construct a full spherical_arc from a plane.
   * \param plane the containing plane.
   */
  Arr_geodesic_arc_on_sphere_3(const Direction_3& normal)
  {
    this->normal(normal);
    this->set_is_vertical(CGAL::sign(normal.dz()) == ZERO);
    this->set_is_directed_right(true);
    this->set_is_full(true);
    this->set_is_degenerate(false);
    this->set_is_empty(false);
    set_is_x_monotone(false);
  }

  /*! Indicates whether the arc is x-monotone
   * \return true if the arc is x-monotone; false otherwise
   */
  bool is_x_monotone() const { return m_is_x_monotone; }

  /*! Set the flag that indicates whether the arc is x-monotone
   * \param flag indicates whether the arc is x-monotone
   */
  void set_is_x_monotone(bool flag) { m_is_x_monotone = flag; }
};

/*! Inserter for the spherical_arc class used by the traits-class */
template <typename Kernel, typename OutputStream>
OutputStream& operator<<(OutputStream& os,
                         const Arr_extended_direction_3<Kernel>& ed)
{
#if defined(CGAL_ARR_GEODESIC_ARC_ON_SPHERE_DETAILS)
  os << "("
     << ed.dx() << ", " << ed.dy() << ",  " << ed.dz();
  os << ")"
     << ", "
     << (ed.is_min_boundary() ? "min" :
         ed.is_max_boundary() ? "max" :
         ed.is_mid_boundary() ? "dis" : "reg");
#else
  // CGAL::To_double<typename Kernel::FT> todouble;
  // os << static_cast<float>(todouble(ed.dx())) << ", "
  //    << static_cast<float>(todouble(ed.dy())) << ", "
  //    << static_cast<float>(todouble(ed.dz()));
#endif
  const typename Kernel::Direction_3* dir = &ed;
  os << *dir;
  return os;
}

/*! Inserter for the spherical_arc class used by the traits-class */
template <typename Kernel, typename OutputStream>
OutputStream&
operator<<(OutputStream& os,
           const Arr_x_monotone_geodesic_arc_on_sphere_3<Kernel>& arc)
{
#if defined(CGAL_ARR_GEODESIC_ARC_ON_SPHERE_DETAILS)
  os << "("
     << "(" << arc.source() << "), (" << arc.target() << ")"
     << "("
     << ", (" << arc.normal() << ")"
     << ", " << (arc.is_vertical() ? " |" : "!|")
     << ", " << (arc.is_directed_right() ? "=>" : "<=")
     << ", " << (arc.is_full() ? "o" : "/");
#else
  os << arc.source() << " " << arc.target() << " ";
  if (!arc.is_meridian()) os << "0";
  else os << "1 " << arc.normal();
#endif
  return os;
}

/*! Extractor for the spherical-arc point class used by the traits-class */
template <typename Kernel, typename InputStream>
InputStream&
operator>>(InputStream& is, Arr_extended_direction_3<Kernel>& point)
{
  typename Kernel::Direction_3* dir = &point;
  is >> *dir;
  point.init();
  return is;
}

/*! Extractor for the spherical_arc class used by the traits-class */
template <typename Kernel, typename InputStream>
InputStream&
operator>>(InputStream& is,
           Arr_x_monotone_geodesic_arc_on_sphere_3<Kernel>& arc)
{
  Arr_extended_direction_3<Kernel> source, target;
  is >> source >> target;
  arc.set_source(source);
  arc.set_target(target);
  unsigned int flag;
  is >> flag;
  if (flag == 1) {
    typename Kernel::Direction_3 normal;
    is >> normal;
    arc.set_normal(normal);
  }
  else {
    Kernel kernel;
    CGAL_precondition(!kernel.equal_3_object()
                      (kernel.construct_opposite_direction_3_object()(source),
                       target));
    typename Kernel::Vector_3 v =
      kernel.construct_cross_product_vector_3_object()(source.vector(),
                                                       target.vector());
    arc.set_normal(v.direction());
  }
  arc.init();
  return is;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>
#endif
