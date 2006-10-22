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

#ifndef CGAL_ARR_BEZIER_CURVE_TRAITS_2_H
#define CGAL_ARR_BEZIER_CURVE_TRAITS_2_H

/*! \file
 * Definition of the Arr_Bezier_curve_traits_2 class. 
 */

#include <CGAL/tags.h>
#include <CGAL/Arr_traits_2/Bezier_curve_2.h>
#include <CGAL/Arr_traits_2/Bezier_point_2.h>
#include <CGAL/Arr_traits_2/Bezier_x_monotone_2.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A traits class for maintaining an arrangement of Bezier curves with
 * rational control points.
 *
 * The class is templated with three parameters: 
 * Rat_kernel A kernel that defines the type of control points.
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices and is used to
 *            represent algebraic numbers.
 * Nt_traits A number-type traits class. This class defines the Rational
 *           number type (should be the same as Rat_kernel::FT) and the
 *           Algebraic number type (should be the same as Alg_kernel::FT)
 *           and supports various operations on them.
 */
template <class Rat_kernel_, class Alg_kernel_, class Nt_traits_>
class Arr_Bezier_curve_traits_2 
{
public:

  typedef Rat_kernel_                            Rat_kernel;
  typedef Alg_kernel_                            Alg_kernel;
  typedef Nt_traits_                             Nt_traits;
  typedef Arr_Bezier_curve_traits_2<Rat_kernel,
                                    Alg_kernel,
                                    Nt_traits>   Self;
 
  typedef typename Nt_traits::Integer            Integer;
  typedef typename Rat_kernel::FT                Rational;
  typedef typename Alg_kernel::FT                Algebraic;

  typedef typename Rat_kernel::Point_2           Rat_point_2;
  typedef typename Alg_kernel::Point_2           Alg_point_2;
  
  // Category tags:
  typedef Tag_true                               Has_left_category;
  typedef Tag_true                               Has_merge_category;
  typedef Tag_false                              Has_infinite_category;

  // Traits-class types:
  typedef _Bezier_curve_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits>             Curve_2;

  typedef _Bezier_x_monotone_2<Rat_kernel,
                               Alg_kernel,
                               Nt_traits>        X_monotone_curve_2;

  typedef _Bezier_point_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits>             Point_2;

private:

  // Type definition for the intersection points mapping.
  typedef typename X_monotone_curve_2::Curve_id           Curve_id;
  typedef typename X_monotone_curve_2::Intersection_point_2
                                                          Intersection_point_2;
  typedef typename X_monotone_curve_2::Intersection_map   Intersection_map;

  Intersection_map  _inter_map;   // Mapping curve pairs to their intersection
                                  // points.

public:

  /*!
   * Default constructor.
   */
  Arr_Bezier_curve_traits_2 ()
  {}

  /// \name Functor definitions.
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
      if (p1.is_same (p2))
        return (EQUAL);

      return (CGAL::compare (p1.x(), p2.x()));
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
      if (p1.is_same (p2))
        return (EQUAL);

      Comparison_result   res = CGAL::compare (p1.x(), p2.x());

      if (res != EQUAL)
        return (res);

      res = CGAL::compare (p1.y(), p2.y());

      if (res == EQUAL)
      {
        // If the two points are the same, merge their originator lists.
        p1.merge_originators (p2);
        p2.merge_originators (p1);
      }
      return (res);
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
      // A rational function can never be vertical:
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
  private:

    Intersection_map&  _inter_map;       // The map of intersection points.

  public:

    /*! Constructor. */
    Compare_y_at_x_2 (const Intersection_map& map) :
      _inter_map (const_cast<Intersection_map&> (map))
    {}

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
      return (cv.point_position (p, _inter_map));
    }
  };

  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return (Compare_y_at_x_2 (_inter_map));
  }

  class Compare_y_at_x_left_2
  {
  private:

    Intersection_map&  _inter_map;       // The map of intersection points.

  public:

    /*! Constructor. */
    Compare_y_at_x_left_2 (const Intersection_map& map) :
      _inter_map (const_cast<Intersection_map&> (map))
    {}

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
      return (cv1.compare_to_left (cv2, p, _inter_map));
    }
  };

  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object () const
  {
    return (Compare_y_at_x_left_2 (_inter_map));
  }

  class Compare_y_at_x_right_2
  {
  private:

    Intersection_map&  _inter_map;       // The map of intersection points.

  public:

    /*! Constructor. */
    Compare_y_at_x_right_2 (const Intersection_map& map) :
      _inter_map (const_cast<Intersection_map&> (map))
    {}

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
      return (cv1.compare_to_right (cv2, p, _inter_map));
    }
  };

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return (Compare_y_at_x_right_2 (_inter_map));
  }

  class Equal_2
  {
  private:

    Intersection_map&  _inter_map;       // The map of intersection points.

  public:

    /*! Constructor. */
    Equal_2 (const Intersection_map& map) :
      _inter_map (const_cast<Intersection_map&> (map))
    {}

    /*!
     * Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      return (cv1.equals (cv2, _inter_map));
    }

    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator() (const Point_2& p1, const Point_2& p2) const
    {
      bool   eq = p1.equals (p2);

      if (eq)
      {
        // If the two points are the same, merge their originator lists.
        p1.merge_originators (p2);
        p2.merge_originators (p1);
      }

      return (eq);
    }
  };

  /*! Get an Equal_2 functor object. */
  Equal_2 equal_2_object () const
  {
    return (Equal_2 (_inter_map));
  }

  class Make_x_monotone_2
  {
  public:

    /*!
     * Cut the given Bezier curve into x-monotone subcurves and insert them
     * into the given output iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The returned
     *           objects is a wrapper for an X_monotone_curve_2 object.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const Curve_2& B, OutputIterator oi)
    {
      // Compute the t-values where B(t) is a point with a vertical tangent.  
      std::list<Algebraic>                           ts;

      B.vertical_tangency_points (std::back_inserter (ts));

      // Create the x-monotone subcurves.
      Algebraic                                      t0 = Algebraic (0);
      typename std::list<Algebraic>::const_iterator  it;

      for (it = ts.begin(); it != ts.end(); ++it)
      {
        *oi = make_object (X_monotone_curve_2 (B, t0, *it));
        ++oi;

        t0 = *it;
      }

      // Create the final subcurve.
      *oi = make_object (X_monotone_curve_2 (B, t0, Algebraic (1)));

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

  class Intersect_2
  {
  private:

    Intersection_map&  _inter_map;       // The map of intersection points.

  public:

    /*! Constructor. */
    Intersect_2 (Intersection_map& map) :
      _inter_map (map)
    {}

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
      return (cv1.intersect (cv2, _inter_map, oi));
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object ()
  {
    return (Intersect_2 (_inter_map));
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
      return (cv1.can_merge_with (cv2));
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
     *      same conic curve and share a common endpoint.
     */
    void operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c) const
    {
      c = cv1.merge (cv2);
      return;
    }
  };

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object () const
  {
    return Merge_2();
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
    Comparison_result operator() (const X_monotone_curve_2& cv)
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
     * Construct an opposite x-monotone curve (with swapped source and target).
     * \param cv The curve.
     * \return The opposite curve.
     */
    X_monotone_curve_2 operator() (const X_monotone_curve_2& cv)
    {
      return (cv.flip());
    }
  };

  /*! Get a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  {
    return Construct_opposite_2();
  }

  //@}
};

CGAL_END_NAMESPACE

#endif
