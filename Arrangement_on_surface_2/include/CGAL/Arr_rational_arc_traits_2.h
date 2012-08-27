// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_RATIONAL_ARC_TRAITS_2_H
#define CGAL_ARR_RATIONAL_ARC_TRAITS_2_H

#define CGAL_DEPRECATED_HEADER "<CGAL/Arr_rational_arc_traits_2.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/Arr_rational_function_traits_2.h>"
#include <CGAL/internal/deprecation_warning.h>

/*! \file
 * Definition of the Arr_rational_arc_traits_2 class. 
 */

#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_geometry_traits/Rational_arc_2.h>

namespace CGAL {

/*! \class
 * A traits class for maintaining an arrangement of bounded arcs (segments) of
 * rational functions of arbitrary degree.
 *
 * The class is templated with two parameters: 
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices, which are algebraic
 *            numbers (defined by Nt_traits::Algebraic).
 * Nt_traits A traits class for performing various operations on the integer,
 *           rational and algebraic types. 
 */
template <class Alg_kernel_, class Nt_traits_>
class Arr_rational_arc_traits_2 
{
public:

  typedef Alg_kernel_                             Alg_kernel;
  typedef Nt_traits_                              Nt_traits;

  typedef typename Alg_kernel::FT                 Algebraic;

  typedef typename Nt_traits::Integer             Integer;

  // Category tags:
  typedef Tag_true                                Has_left_category;
  typedef Tag_true                                Has_merge_category;
  typedef Tag_false                               Has_do_intersect_category;

  typedef Arr_open_side_tag                       Left_side_category;
  typedef Arr_open_side_tag                       Bottom_side_category;
  typedef Arr_open_side_tag                       Top_side_category;
  typedef Arr_open_side_tag                       Right_side_category;

  // Traits objects:
  typedef _Rational_arc_2<Alg_kernel, Nt_traits>  Curve_2;
  typedef _Continuous_rational_arc_2<Alg_kernel, Nt_traits>
                                                  X_monotone_curve_2;
  typedef typename Alg_kernel::Point_2            Point_2;
  typedef unsigned int                            Multiplicity;

  typedef typename Curve_2::Rat_vector            Rat_vector;

public:

  /*!
   * Default constructor.
   */
  Arr_rational_arc_traits_2 ()
  {}

  /// \name Functor definitions.
  //@{

  /*! A functor that compares the x-coordinates of two points */
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

  /*! Obtain a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object () const
  {
    return Compare_x_2();
  }

  /*! A functor that compares two points lexigoraphically: by x, then by y. */
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

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return Compare_xy_2();
  }

  /*! A functor that obtains the left endpoint of a curve. */
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

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object () const
  {
    return Construct_min_vertex_2();
  }

  /*! A functor that obtains the right endpoint of a curve. */
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

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object () const
  {
    return Construct_max_vertex_2();
  }

  /*! A functor that checks whether a given curve is vertical. */
  class Is_vertical_2
  {
  public:
    /*!
     * Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2&) const
    {
      // A rational function can never be vertical.
      return (false);
    }
  };

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object () const
  {
    return Is_vertical_2();
  }

  /*! A functor that compares the y-coordinates of a point and a curve at
   * the point x-coordinate.
   */
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
      return (cv.point_position (p));
    }
  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return Compare_y_at_x_2();
  }

  /*! A functor that compares compares the y-coordinates of two curves
   * immediately to the left of their intersection point.
   */
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
      CGAL_precondition (cv1.point_position (p) == EQUAL &&
                         cv2.point_position (p) == EQUAL);

      CGAL_precondition_code (
        Alg_kernel   ker;
      );
      CGAL_precondition ((cv1.left_infinite_in_x() != ARR_INTERIOR ||
                          cv1.left_infinite_in_y() != ARR_INTERIOR ||
                          ker.compare_xy_2_object() (p, 
                                                     cv1.left()) == LARGER) &&
                         (cv2.left_infinite_in_x() != ARR_INTERIOR ||
                          cv2.left_infinite_in_y() != ARR_INTERIOR ||
                          ker.compare_xy_2_object() (p,
                                                     cv2.left()) == LARGER));

      // Compare the slopes of the two arcs.
      Comparison_result        res;
      unsigned int             mult;
    
      res = cv1.compare_slopes (cv2, p, mult);

      // The comparison result is to the right of p. In case the multiplicity 
      // of the intersection point p is odd, reverse this result.
      if (mult % 2 == 1)
      {
	if (res == SMALLER)
	  res = LARGER;
	else if (res == LARGER)
	  res = SMALLER;
      }

      return (res);
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object () const
  {
    return Compare_y_at_x_left_2();
  }

  /*! A functor that compares compares the y-coordinates of two curves
   * immediately to the right of their intersection point.
   */
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
      CGAL_precondition (cv1.point_position (p) == EQUAL &&
                         cv2.point_position (p) == EQUAL);

      CGAL_precondition_code (
        Alg_kernel   ker;
      );
      CGAL_precondition((cv1.right_infinite_in_x() != ARR_INTERIOR ||
                         cv1.right_infinite_in_y() != ARR_INTERIOR ||
                         ker.compare_xy_2_object() (p, 
                                                    cv1.right()) == SMALLER) &&
                        (cv2.right_infinite_in_x() != ARR_INTERIOR ||
                         cv2.right_infinite_in_y() != ARR_INTERIOR ||
                         ker.compare_xy_2_object() (p,
                                                    cv2.right()) == SMALLER));

      // Compare the slopes of the two arcs to determine thir relative
      // position immediately to the right of p.
      unsigned int             mult;

      return (cv1.compare_slopes (cv2, p, mult));
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return Compare_y_at_x_right_2();
  }

  /*! A functor that checks whether two points and two curves are identical. */
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

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object () const
  {
    return Equal_2();
  }

  /*! A functor that divides a curve into continues (x-monotone) curves. */
  class Make_x_monotone_2
  {
  public:

    /*!
     * Cut the given conic curve (or conic arc) into x-monotone subcurves 
     * and insert them to the given output iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The returned
     *           objects is a wrapper for an X_monotone_curve_2 object.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi)
    {
      // Make the rational arc continuous.
      std::list<X_monotone_curve_2>                           arcs;

      cv.make_continuous (std::back_inserter (arcs));

      // Create objects.
      typename std::list<X_monotone_curve_2>::const_iterator  iter;

      for (iter = arcs.begin(); iter != arcs.end(); ++iter)
      {
        *oi = make_object (*iter);
        ++oi;
      }

      return (oi);
    }
  };

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () const
  {
    return Make_x_monotone_2();
  }

  /*! A functor that splits a curve at a point. */
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

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object () const
  {
    return Split_2();
  }

  /*! A functor that computes intersections between two curves. */
  class Intersect_2
  {
  public:

    /*!
     * Find the intersections of the two given curves and insert them to the
     * given output iterator. As two segments may itersect only once, only a
     * single will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const X_monotone_curve_2& cv1,
                               const X_monotone_curve_2& cv2,
                               OutputIterator oi)
    {
      return (cv1.intersect (cv2, oi));
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return Intersect_2();
  }

  /*! A functor that tests whether two curves can be merged. */
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

  /*! Obtain an Are_mergeable_2 functor object. */
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
    typedef Arr_rational_arc_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_rational_arc_traits_2<Kernel>;
    
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
      CGAL_precondition(m_traits->are_mergeable_2_object()(cv1, cv2));

      c = cv1;
      c.merge (cv2);
    }
  };

  /*! Obtain a Merge_2 functor object. */
  Merge_2 merge_2_object () const
  {
    return Merge_2(this);
  }

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
      return (ce == ARR_MIN_END) ?
        xcv.left_infinite_in_x() : xcv.right_infinite_in_x();
    }

    /*! Obtains the parameter space at a point along the x-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    Arr_parameter_space operator()(const Point_2 ) const
    {
      return ARR_INTERIOR;
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
      return (ce == ARR_MIN_END) ?
        xcv.left_infinite_in_y() : xcv.right_infinite_in_y();
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     */
    Arr_parameter_space operator()(const Point_2 ) const
    {
      return ARR_INTERIOR;
    }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(); }


 /*! A function object that compares the x-limits of arc ends on the
   * boundary of the parameter space
   */
  class Compare_x_at_limit_2 {
  public:
    /*! Compare the x-coordinate of a point and the x-coordinate of the limit
     * of a rational arc at its specificed end at y = +/- oo.
     * \param p the point.
     * \param xcv the parametric curve, the endpoint of which is compared.
     * \param ce identifies an open end of xcv -
     *            ARR_MIN_END - the minimal end of xc or
     *            ARR_MAX_END - the maximal end of xc.
     * \return the comparison result:
     *         SMALLER - x(p) < x(xc, ce);
     *         EQUAL   - x(p) = x(xc, ce);
     *         LARGER  - x(p) > x(xc, ce).     
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the curve xcv lies on a boundary, implying that xcv
     *      is either a vertical line or a curve with a vertical asymptote.
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & xcv,
                                 Arr_curve_end ce) const
    {
      // TODO implement (simplify)
      return (xcv.compare_end (p, ce));
    }

    /*! Compare the x-coordinates of the limits of 2 parameteric curves at
     * their specificed ends at y = +/- oo.
     * \param xcv1 the first arc.
     * \param ce1 identifies an open end of ce1 -
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second arc.
     * \param ce2 identifies an open end of ce2 -
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce1 end of the line xcv1 lies on a boundary, implying that xcv1
     *      is either a vertical line or a curve with a vertical asymptote.
     * \pre the ce2 end of the line xcv2 lies on a boundary, implying that xcv2
     *      is either a vertical line or a curve with a vertical asymptote.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce2) const
    {
      // TODO implement (simplify)
      return xcv1.compare_ends (ce1, xcv2, ce2);
    }
  };

  /*! Obtain a Compare_x_at_limit_2 function object */
  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(); }
    
  /*! A function object that compares the x-coordinates of arc ends near the
   * boundary of the parameter space
   */
  class Compare_x_near_limit_2 {
  public:
    
    /*! Compare the x-coordinates of 2 arcs ends near the boundary of the
     * parameter space at y = +/- oo.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the curve end indicator -
     *            ARR_MIN_END - the minimal end of curves or
     *            ARR_MAX_END - the maximal end of curves.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce) < x(xcv2, ce);
     *         EQUAL   - x(xcv1, ce) = x(xcv2, ce);
     *         LARGER  - x(xcv1, ce) > x(xcv2, ce).
     * \pre the ce end of the line xcv1 lies on a boundary.
     * \pre the ce end of the line xcv2 lies on a boundary.
     */
    Comparison_result operator()(const X_monotone_curve_2 & xcv1,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce) const
    {
      // TODO implement (simplify)
      return xcv1.compare_ends (ce, xcv2, ce);
    }
  };

  /*! Obtain a Compare_x_near_limit_2 function object */
  Compare_x_near_limit_2 compare_x_near_limit_2_object() const
  { return Compare_x_near_limit_2(); }
    

  /*! A function object that compares the y-coordinates of arc ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  public:
    /*! Compare the y-coordinates of 2 lines at their ends near the boundary
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
      return (ce == ARR_MIN_END) ?
        xcv1.compare_at_minus_infinity (xcv2) :
        xcv1.compare_at_plus_infinity (xcv2);
    }
  };

  /*! Obtain a Compare_y_near_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(); }

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
    Comparison_result operator() (const X_monotone_curve_2& cv)
    {
      if (cv.is_directed_right())
        return (SMALLER);
      else
	return (LARGER);
    }
  };

  /*! Obtain a Compare_endpoints_xy_2 functor object. */
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
    X_monotone_curve_2 operator() (const X_monotone_curve_2& cv)
    {
      return (cv.flip());
    }
  };

  /*! Obtain a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  {
    return Construct_opposite_2();
  }

  //@}
};

} //namespace CGAL

#endif
