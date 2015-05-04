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
// Author(s)     : Ron Wein     <wein@post.tau.ac.il>
//                 Iddo Hanniel <iddoh@cs.technion.ac.il>
//                 Waqar Khan   <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_ARR_BEZIER_CURVE_TRAITS_2_H
#define CGAL_ARR_BEZIER_CURVE_TRAITS_2_H

/*! \file
 * Definition of the Arr_Bezier_curve_traits_2 class.
 */

#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_geometry_traits/Bezier_curve_2.h>
#include <CGAL/Arr_geometry_traits/Bezier_point_2.h>
#include <CGAL/Arr_geometry_traits/Bezier_x_monotone_2.h>
#include <CGAL/Arr_geometry_traits/Bezier_cache.h>
#include <CGAL/Arr_geometry_traits/Bezier_bounding_rational_traits.h>

namespace CGAL {

/*! \class
 * A traits class for maintaining an arrangement of Bezier curves with
 * rational control points.
 *
 * The class is templated with four parameters:
 * Rat_kernel A kernel that defines the type of control points.
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices and is used to
 *            represent algebraic numbers.
 * Nt_traits A number-type traits class. This class defines the Rational
 *           number type (should be the same as Rat_kernel::FT) and the
 *           Algebraic number type (should be the same as Alg_kernel::FT)
 *           and supports various operations on them.
 * Bounding_traits A traits class for filtering the exact computations.
 *                 By default we use the rational bounding traits.
 */
template <class RatKernel_, class AlgKernel_, class NtTraits_,
          class BoundingTraits_ = Bezier_bounding_rational_traits<RatKernel_> >
class Arr_Bezier_curve_traits_2
{
public:

  typedef RatKernel_                             Rat_kernel;
  typedef AlgKernel_                             Alg_kernel;
  typedef NtTraits_                              Nt_traits;
  typedef BoundingTraits_                        Bounding_traits;
  typedef Arr_Bezier_curve_traits_2<Rat_kernel,
                                    Alg_kernel,
                                    Nt_traits,
                                    Bounding_traits>   Self;

  typedef typename Nt_traits::Integer            Integer;
  typedef typename Rat_kernel::FT                Rational;
  typedef typename Alg_kernel::FT                Algebraic;

  typedef typename Rat_kernel::Point_2           Rat_point_2;
  typedef typename Alg_kernel::Point_2           Alg_point_2;

  // Category tags:
  typedef Tag_true                               Has_left_category;
  typedef Tag_true                               Has_merge_category;
  typedef Tag_false                              Has_do_intersect_category;

  typedef Arr_oblivious_side_tag                 Left_side_category;
  typedef Arr_oblivious_side_tag                 Bottom_side_category;
  typedef Arr_oblivious_side_tag                 Top_side_category;
  typedef Arr_oblivious_side_tag                 Right_side_category;

  // Traits-class types:
  typedef _Bezier_curve_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits,
                          Bounding_traits>             Curve_2;

  typedef _Bezier_x_monotone_2<Rat_kernel,
                               Alg_kernel,
                               Nt_traits,
                               Bounding_traits>        X_monotone_curve_2;

  typedef _Bezier_point_2<Rat_kernel,
                          Alg_kernel,
                          Nt_traits,
                          Bounding_traits>             Point_2;

  typedef typename X_monotone_curve_2::Multiplicity    Multiplicity;


  // Type definition for the vertical-tangnecy and intersection point cache.
  typedef _Bezier_cache<Nt_traits>                     Bezier_cache;

private:

  // Type definition for the bounded intersection points mapping.
  typedef typename X_monotone_curve_2::Intersection_map   Intersection_map;

  // Data members:
  mutable Bezier_cache * p_cache;         /*!< Caches vertical tangency points
                                           * and intersection points that have
                                           * been computed in an exact manner.
                                           */
  mutable Intersection_map * p_inter_map; /*!< Maps curve pairs to their
                                           * intersection points.
                                           */
  bool m_owner;                           /*!< Does this instance own its cache
                                           * and map structures.
                                           */

public:

  /// \name Construction.
  //@{

  /*! Default constructor. */
  Arr_Bezier_curve_traits_2 ()
  {
    p_cache = new Bezier_cache;
    p_inter_map = new Intersection_map;
    m_owner = true;
  }

  /*! Copy constructor. */
  Arr_Bezier_curve_traits_2 (const Self& tr) :
    p_cache (tr.p_cache),
    p_inter_map (tr.p_inter_map),
    m_owner (false)
  {}

  /*! Assignmnet operator. */
  Self& operator= (const Self& tr)
  {
    if (this == &tr)
      return (*this);

    p_cache = tr.p_cache;
    p_inter_map = tr.p_inter_map;
    m_owner = false;
    return (*this);
  }

  /*! Destructor. */
  ~Arr_Bezier_curve_traits_2 ()
  {
    if (m_owner)
    {
      delete p_cache;
      delete p_inter_map;
    }
    p_cache = NULL;
    p_inter_map = NULL;
  }
  //@}

  /// \name Functor definitions for the arrangement traits.
  //@{

  /*! \class Compare_x_2
   * The Compare_x_2 functor.
   */
  class Compare_x_2
  {
  private:
    const Bezier_cache   *p_cache;

  public:

    /*! Constructor. */
    Compare_x_2 (const Bezier_cache *cache) :
      p_cache (cache)
    {}

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
      return (p1.compare_x (p2,
                            const_cast<Bezier_cache&> (*p_cache)));
    }
  };

  /*! Get a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object () const
  {
    return (Compare_x_2 (p_cache));
  }

  /*! \class Compare_xy_2
   * The Compare_xy_2 functor.
   */
  class Compare_xy_2
  {
  private:
    const Bezier_cache   *p_cache;

  public:

    /*! Constructor. */
    Compare_xy_2 (const Bezier_cache *cache) :
      p_cache (cache)
    {}

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
      return (p1.compare_xy (p2,
                             const_cast<Bezier_cache&> (*p_cache)));
    }
  };

  /*! Get a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object () const
  {
    return (Compare_xy_2 (p_cache));
  }

  /*! \class Construct_min_vertex_2
   * The Construct_min_vertex_2 functor.
   */
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

  /*! \class Construct_max_vertex_2
   * The Construct_max_vertex_2 functor.
   */
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

  /*! \class Is_vertical_2
   * The Is_vertical_2 functor.
   */
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

  /*! \class Compare_y_at_x_2
   * The Compare_y_at_x_2 functor.
   */
  class Compare_y_at_x_2
  {
  private:
    const Bezier_cache   *p_cache;

  public:

    /*! Constructor. */
    Compare_y_at_x_2 (const Bezier_cache *cache) :
      p_cache (cache)
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
      return (cv.point_position (p,
                                 const_cast<Bezier_cache&> (*p_cache)));
    }
  };

  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return (Compare_y_at_x_2 (p_cache));
  }

  /*! \class Compare_y_at_x_left_2
   * The Compare_y_at_x_left_2 functor.
   */
  class Compare_y_at_x_left_2
  {
  private:
    const Bezier_cache   *p_cache;

  public:

    /*! Constructor. */
    Compare_y_at_x_left_2 (const Bezier_cache *cache) :
      p_cache (cache)
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
      return (cv1.compare_to_left (cv2, p,
                                   const_cast<Bezier_cache&> (*p_cache)));
    }
  };

  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object () const
  {
    return (Compare_y_at_x_left_2 (p_cache));
  }

  /*! \class Compare_y_at_x_right_2
   * The Compare_y_at_x_right_2 functor.
   */
  class Compare_y_at_x_right_2
  {
  private:
    const Bezier_cache   *p_cache;

  public:

    /*! Constructor. */
    Compare_y_at_x_right_2 (const Bezier_cache *cache) :
      p_cache (cache)
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
      return (cv1.compare_to_right (cv2, p,
                                    const_cast<Bezier_cache&> (*p_cache)));
    }
  };

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object () const
  {
    return (Compare_y_at_x_right_2 (p_cache));
  }

  /*! \class Equal_2
   * The Equal_2 functor.
   */
  class Equal_2
  {
  private:
    const Bezier_cache         *p_cache;

  public:

    /*! Constructor. */
    Equal_2 (const Bezier_cache *cache) :
      p_cache (cache)
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
      return (cv1.equals (cv2,
                          const_cast<Bezier_cache&> (*p_cache)));
    }

    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    bool operator() (const Point_2& p1, const Point_2& p2) const
    {
      return (p1.equals (p2,
                         const_cast<Bezier_cache&> (*p_cache)));
    }
  };

  /*! Get an Equal_2 functor object. */
  Equal_2 equal_2_object () const
  {
    return (Equal_2 (p_cache));
  }

  /*! \class Make_x_monotone_2
   * The Make_x_monotone_2 functor.
   */
  class Make_x_monotone_2
  {
  private:
    Bezier_cache         *p_cache;

  public:

    /*! Constructor. */
    Make_x_monotone_2 (Bezier_cache *cache) :
      p_cache (cache)
    {}

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
      typedef typename Bounding_traits::Vertical_tangency_point
                                                Vertical_tangency_point;

      // Try to compute the bounds of the vertical tangency points.
      Bounding_traits                           bound_tr;
      typename Bounding_traits::Control_points  cpts;
      std::list<Vertical_tangency_point>        vpt_bounds;

      std::copy (B.control_points_begin(), B.control_points_end(),
                 std::back_inserter(cpts));

      bound_tr.compute_vertical_tangency_points
          (cpts, std::back_inserter (vpt_bounds));

      // Construct Point_2 from bounded tangency points.
      std::list<Point_2>                            vpts;
      bool                                          app_ok = true;
      typename std::list<Vertical_tangency_point>::const_iterator iter;

      for (iter = vpt_bounds.begin(); iter != vpt_bounds.end(); ++iter)
      {
        const typename Bounding_traits::Bez_point_bound& bound = iter->bound;
        const typename Bounding_traits::Bez_point_bbox&  bbox = iter->bbox;

        if (! bound.can_refine)
        {
          // If we cannot refine the vertical-tangency bound anymore, then
          // we failed to bound the vertical tangency point.
          // \todo In the future, we might want to use the info.
          app_ok = false;
          break;
        }

        // Construct an approximate vertical tangency point.
        Point_2   pt;

        if (bound.type == Bounding_traits::Bez_point_bound::RATIONAL_PT)
        {
          CGAL_assertion (CGAL::compare (bound.t_min, bound.t_max) == EQUAL);
          Rational  t0 = bound.t_min;

          pt = Point_2 (B, t0);
        }
        else
        {
          pt.add_originator (typename Point_2::Originator(B, bound));
        }
        pt.set_bbox (bbox);

        vpts.push_back(pt);
      }

      // If bounding the approximated vertical-tangency points went fine,
      // use these points as endpoint for the x-monotone subcurves.
      if (app_ok)
      {
        // Create the x-monotone subcurves with approximate endpoints.
        typename std::list<Point_2>::const_iterator pit;
        Point_2       p0(B, Rational(0)); // A rational start point.
        unsigned int  xid = 1;            // Serial number of the subcurve.

        for (pit = vpts.begin(); pit != vpts.end(); ++pit)
        {
          *oi++ = CGAL::make_object (X_monotone_curve_2 (B, xid,
                                                         p0, *pit,
                                                         *p_cache));
          xid++;
          p0 = *pit;
        }

        Point_2    p1(B, Rational(1)); // A rational end point.

        *oi++ = CGAL::make_object (X_monotone_curve_2 (B, xid,
                                                       p0, p1,
                                                       *p_cache));
        return (oi);
      }

      // If we reached here then we have to compute the vertical-tangency
      // points in an exact manner. We do this by considering all t-values
      // on B(t) = (X(t), Y(t)), such that X'(t) = 0.
      const typename Bezier_cache::Vertical_tangency_list&
        vt_list = p_cache->get_vertical_tangencies (B.id(),
                                                    B.x_polynomial(),
                                                    B.x_norm());

      // Create the x-monotone subcurves.
      Point_2                                        p0 (B, Rational(0));
      Point_2                                        p1;
      typename Bezier_cache::Vertical_tangency_iter  it;
      unsigned int  xid = 1;            // Serial number of the subcurve.

      for (it = vt_list.begin(); it != vt_list.end(); ++it)
      {
        p1 = Point_2 (B, *it);
        *oi++ = CGAL::make_object (X_monotone_curve_2 (B, xid,
                                                       p0, p1,
                                                       *p_cache));
        xid++;
        p0 = p1;
      }

      // Create the final subcurve.
      p1 = Point_2 (B, Rational(1));
      *oi++ = CGAL::make_object (X_monotone_curve_2 (B, xid,
                                                     p0, p1,
                                                     *p_cache));
      return (oi);
    }
  };

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () const
  {
    return (Make_x_monotone_2 (p_cache));
  }

  /*! \class Split_2
   * The Split_2 functor.
   */
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

  /*! \class Intersect_2
   * The Intersect_2 functor.
   */
  class Intersect_2
  {
  private:
    Bezier_cache         *p_cache;
    Intersection_map     *p_imap;

  public:

    /*! Constructor. */
    Intersect_2 (Bezier_cache *cache, Intersection_map *imap) :
      p_cache (cache),
      p_imap (imap)
    {}

    /*!
     * Find the intersections of the two given curves and insert them to the
     * given output iterator.
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
      return (cv1.intersect (cv2, *p_imap, *p_cache, oi));
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return (Intersect_2 (p_cache, p_inter_map));
  }

  /*! \class Are_mergeable_2
   * The Are_mergeable_2 functor.
   */
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
    typedef Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel,
                                      Nt_traits, Bounding_traits>       Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel,
                                           Nt_traits, Bounding_traits>;

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

      c = cv1.merge (cv2);
      return;
    }
  };

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object () const
  {
    return Merge_2(this);
  }

  //@}

  /// \name Functor definitions for the Boolean set-operation traits.
  //@{

  /*! \class Compare_endpoints_xy_2
   * The Compare_endpoints_xy_2 functor.
   */
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

  class Trim_2 {
    typedef Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel,
                                      Nt_traits, Bounding_traits>       Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Trim_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel,
                                           Nt_traits, Bounding_traits>;
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
  public:
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv,
                                  const Point_2& src,
                                  const Point_2& tgt) const
    {
      // make functor objects
      CGAL_precondition_code(Compare_y_at_x_2 compare_y_at_x_2 =
                             m_traits.compare_y_at_x_2_object());
      CGAL_precondition_code(Equal_2 equal_2 = m_traits.equal_2_object());
      Compare_x_2 compare_x_2 = m_traits.compare_x_2_object();
      // Check whether source and taget are two distinct points and they lie
      // on the line.
      CGAL_precondition(compare_y_at_x_2(src, xcv) == EQUAL);
      CGAL_precondition(compare_y_at_x_2(tgt, xcv) == EQUAL);
      CGAL_precondition(! equal_2(src, tgt));

      //check if the orientation conforms to the src and tgt.
      if( xcv.is_directed_right() && compare_x_2(src, tgt) == LARGER)
        return (xcv.trim(tgt, src));
      else if(! xcv.is_directed_right() && compare_x_2(src, tgt) == SMALLER)
        return (xcv.trim(tgt, src));
      else return (xcv.trim(src, tgt));
    }
  };

  /*! Obtain a Trim_2 functor object. */
  Trim_2 trim_2_object() const { return Trim_2(*this); }

  /*! \class Construct_opposite_2
   * The Construct_opposite_2 functor.
   */
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

} //namespace CGAL

#endif
