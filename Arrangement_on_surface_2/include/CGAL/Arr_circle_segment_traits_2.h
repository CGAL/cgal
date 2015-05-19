// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Baruch Zukerman   <baruchzu@post.tau.ac.il>
//                 Waqar Khan <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_ARR_CIRCLE_SEGMENT_TRAITS_2_H
#define CGAL_ARR_CIRCLE_SEGMENT_TRAITS_2_H

/*! \file
 * The header file for the Arr_circle_segment_traits_2<Kenrel> class.
 */

#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_geometry_traits/Circle_segment_2.h>

#include <fstream>

namespace CGAL {

/*! \class
 * A traits class for maintaining an arrangement of circles.
 */
template <typename Kernel_, bool Filter = true>
class Arr_circle_segment_traits_2 {
public:
  typedef Kernel_                                        Kernel;
  typedef typename Kernel::FT                            NT;
  typedef typename Kernel::Point_2                       Rational_point_2;
  typedef typename Kernel::Segment_2                     Rational_segment_2;
  typedef typename Kernel::Circle_2                      Rational_circle_2;
  typedef _One_root_point_2<NT, Filter>                  Point_2;
  typedef typename Point_2::CoordNT                      CoordNT;
  typedef _Circle_segment_2<Kernel, Filter>              Curve_2;
  typedef _X_monotone_circle_segment_2<Kernel, Filter>   X_monotone_curve_2;
  typedef unsigned int                                   Multiplicity;
  typedef Arr_circle_segment_traits_2<Kernel, Filter>    Self;

  // Category tags:
  typedef Tag_true                                   Has_left_category;
  typedef Tag_true                                   Has_merge_category;
  typedef Tag_false                                  Has_do_intersect_category;

  typedef Arr_oblivious_side_tag                     Left_side_category;
  typedef Arr_oblivious_side_tag                     Bottom_side_category;
  typedef Arr_oblivious_side_tag                     Top_side_category;
  typedef Arr_oblivious_side_tag                     Right_side_category;

protected:
  // Type definition for the intersection points mapping.
  typedef typename X_monotone_curve_2::Intersection_map   Intersection_map;

  mutable Intersection_map inter_map;   // Mapping pairs of curve IDs to their
                                        // intersection points.
  bool m_use_cache;

public:
  /*! Default constructor. */
  Arr_circle_segment_traits_2 (bool use_intersection_caching = false) :
    m_use_cache(use_intersection_caching)
  {}

  /*! Get the next curve index. */
  static unsigned int get_index ()
  {
    static unsigned int index = 0;
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
    Comparison_result operator() (const Point_2& p1, const Point_2& p2) const
    {
      if (p1.identical (p2))
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
      if (p1.identical (p2))
        return (EQUAL);

      Comparison_result  res = CGAL::compare (p1.x(), p2.x());

      if (res != EQUAL)
        return (res);

      return (CGAL::compare (p1.y(), p2.y()));
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
    Comparison_result operator() (const Point_2& p,
                                  const X_monotone_curve_2& cv) const
    {
      CGAL_precondition (cv.is_in_x_range (p));

      return (cv.point_position (p));
    }
  };

  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object () const
  {
    return Compare_y_at_x_2();
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
      // right (so their right endpoint is lexicographically larger than p).
      CGAL_precondition (cv1.point_position (p) == EQUAL &&
                         cv2.point_position (p) == EQUAL);

      if ((CGAL::compare (cv1.left().x(),cv1.right().x()) == EQUAL) &&
          (CGAL::compare (cv2.left().x(),cv2.right().x()) == EQUAL))
      { //both cv1 and cv2 are vertical
        CGAL_precondition (!(cv1.right()).equals(p) && !(cv2.right()).equals(p));
      }
      else if ((CGAL::compare (cv1.left().x(),cv1.right().x()) != EQUAL) &&
               (CGAL::compare (cv2.left().x(),cv2.right().x()) == EQUAL))
      { //only cv1 is vertical
        CGAL_precondition (!(cv1.right()).equals(p));
      }
      else if ((CGAL::compare (cv1.left().x(),cv1.right().x()) == EQUAL) &&
               (CGAL::compare (cv2.left().x(),cv2.right().x()) != EQUAL))
      { //only cv2 is vertical
        CGAL_precondition (!(cv2.right()).equals(p));
      }
      else
      { //both cv1 and cv2 are non vertical
        CGAL_precondition (CGAL::compare (cv1.right().x(),p.x()) == LARGER &&
                           CGAL::compare (cv2.right().x(),p.x()) == LARGER);
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

      if ((CGAL::compare (cv1.left().x(),cv1.right().x()) == EQUAL) &&
          (CGAL::compare (cv2.left().x(),cv2.right().x()) == EQUAL))
	  { //both cv1 and cv2 are vertical
         CGAL_precondition (!(cv1.left()).equals(p) && !(cv2.left()).equals(p));
	  }
	  else if ((CGAL::compare (cv1.left().x(),cv1.right().x()) != EQUAL) &&
                   (CGAL::compare (cv2.left().x(),cv2.right().x()) == EQUAL))
	  { //only cv1 is vertical
         CGAL_precondition (!(cv1.left()).equals(p));
	  }
	  else if ((CGAL::compare (cv1.left().x(),cv1.right().x()) == EQUAL) &&
                   (CGAL::compare (cv2.left().x(),cv2.right().x()) != EQUAL))
	  { //only cv2 is vertical
         CGAL_precondition (!(cv2.left()).equals(p));
	  }
	  else
	  { //both cv1 and cv2 are non vertical
        CGAL_precondition (CGAL::compare (cv1.left().x(),p.x()) == SMALLER &&
                           CGAL::compare (cv2.left().x(),p.x()) == SMALLER);
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
      return (p1.equals (p2));
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
  private:
    typedef Arr_circle_segment_traits_2<Kernel_, Filter> Self;

    bool m_use_cache;

  public:

    Make_x_monotone_2(bool use_cache = false) : m_use_cache(use_cache)
    {}

    /*!
     * Cut the given conic curve (ocv.is_in_x_range (p)r conic arc) into x-monotone subcurves
     * and insert them to the given output iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The returned
     *           objects are all wrcv.is_in_x_range (p)appers X_monotone_curve_2 objects.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi) const
    {
      // Increment the serial number of the curve cv, which will serve as its
      // unique identifier.
      unsigned int  index = 0;
      if(m_use_cache)
        index = Self::get_index();

      if (cv.orientation() == COLLINEAR)
      {
        // The curve is a line segment.
        *oi = make_object (X_monotone_curve_2 (cv.supporting_line(),
                                               cv.source(), cv.target(),
                                               index));
        ++oi;
        return (oi);
      }

      // Check the case of a degenrate circle (a point).
      const typename Kernel::Circle_2&  circ = cv.supporting_circle();
      CGAL::Sign   sign_rad = CGAL::sign (circ.squared_radius());
      CGAL_precondition (sign_rad != NEGATIVE);

      if (sign_rad == ZERO)
      {
        // Create an isolated point.
        *oi = make_object (Point_2 (circ.center().x(), circ.center().y()));
        ++oi;
        return (oi);
      }

      // The curve is circular: compute the to vertical tangency points
      // of the supporting circle.
      Point_2         vpts[2];
      unsigned int    n_vpts = cv.vertical_tangency_points (vpts);

      if (cv.is_full())
      {
        CGAL_assertion (n_vpts == 2);

        // Subdivide the circle into two arcs (an upper and a lower half).
        *oi = make_object (X_monotone_curve_2 (circ,
                                               vpts[0], vpts[1],
                                               cv.orientation(),
                                               index));
        ++oi;

        *oi = make_object (X_monotone_curve_2 (circ,
                                               vpts[1], vpts[0],
                                               cv.orientation(),
                                               index));
        ++oi;
      }
      else
      {
        // Act according to the number of vertical tangency points.
        if (n_vpts == 2)
        {
          // Subdivide the circular arc into three x-monotone arcs.
          *oi = make_object (X_monotone_curve_2 (circ,
                                                 cv.source(), vpts[0],
                                                 cv.orientation(),
                                                 index));
          ++oi;

          *oi = make_object (X_monotone_curve_2 (circ,
                                                 vpts[0], vpts[1],
                                                 cv.orientation(),
                                                 index));
          ++oi;

          *oi = make_object (X_monotone_curve_2 (circ,
                                                 vpts[1], cv.target(),
                                                 cv.orientation(),
                                                 index));
          ++oi;
        }
        else if (n_vpts == 1)
        {
          // Subdivide the circular arc into two x-monotone arcs.
          *oi = make_object (X_monotone_curve_2 (circ,
                                                 cv.source(), vpts[0],
                                                 cv.orientation(),
                                                 index));
          ++oi;

          *oi = make_object (X_monotone_curve_2 (circ,
                                                 vpts[0], cv.target(),
                                                 cv.orientation(),
                                                 index));
          ++oi;
        }
        else
        {
          CGAL_assertion (n_vpts == 0);

          // The arc is already x-monotone:
          *oi = make_object (X_monotone_curve_2 (circ,
                                                 cv.source(), cv.target(),
                                                 cv.orientation(),
                                                 index));
          ++oi;
        }
      }

      return (oi);
    }
  };

  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object () const
  {
    return Make_x_monotone_2(m_use_cache);
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
      CGAL_precondition (cv.point_position(p)==EQUAL &&
      ! p.equals (cv.source()) &&
      ! p.equals (cv.target()));

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
                               OutputIterator oi) const
    {
      return (cv1.intersect (cv2, oi, &_inter_map));
    }
  };

  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object () const
  {
    return (Intersect_2 (inter_map));
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

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2
  {
  protected:
    typedef Arr_circle_segment_traits_2<Kernel, Filter> Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_circle_segment_traits_2<Kernel, Filter>;

  public:
    /*!
     * Merge two given x-monotone curves into a single curve.
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

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object () const
  {
    return Merge_2(this);
  }

  class Compare_endpoints_xy_2
  {
  public:
    /*!
     * compare lexicogrphic the endpoints of a x-monotone curve.
     * \param cv the curve
     * \return SMALLER if the curve is directed right, else return SMALLER
     */
    Comparison_result operator()(const X_monotone_curve_2& cv) const
    {
      if(cv.is_directed_right())
        return(SMALLER);
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
     * construct an opposite x-monotone curve.
     * \param cv the curve
     * \return an opposite x-monotone curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& cv) const
    {
      return cv.construct_opposite();
    }
  };

  /*! Get a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  {
    return Construct_opposite_2();
  }

  class Trim_2 {
  protected:
    typedef Arr_circle_segment_traits_2<Kernel, Filter> Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Trim_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_circle_segment_traits_2<Kernel, Filter>;

  public:
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
      // Check whether source and taget are two distinct points and they lie
      // on the line.
      CGAL_precondition(compare_y_at_x_2(src, xcv) == EQUAL);
      CGAL_precondition(compare_y_at_x_2(tgt, xcv) == EQUAL);
      CGAL_precondition(! equal_2(src, tgt));

      //check if the orientation conforms to the src and tgt.
      if( (xcv.is_directed_right() && compare_x_2(src, tgt) == LARGER) ||
          (! xcv.is_directed_right() && compare_x_2(src, tgt) == SMALLER) )
        return (xcv.trim(tgt, src) );
      else return (xcv.trim(src, tgt));
    }
  };

  /*! Obtain a Trim_2 functor object. */
  Trim_2 trim_2_object() const { return Trim_2(*this); }
};

} //namespace CGAL

#endif
