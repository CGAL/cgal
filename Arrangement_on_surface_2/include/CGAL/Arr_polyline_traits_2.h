// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University(Israel).
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
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>
//                 Ron Wein  <wein@post.tau.ac.il>

#ifndef CGAL_ARR_POLYLINE_TRAITS_2_H
#define CGAL_ARR_POLYLINE_TRAITS_2_H

/*! \file
 * The traits-class for the linear piece-wiese(polyline) type of curves of the
 * arrangement package.
 */

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_geometry_traits/Polyline_2.h>

namespace CGAL {

template <class T_SegmentTraits_2>
class Arr_polyline_traits_2 {
public:
  typedef T_SegmentTraits_2                          Segment_traits_2;

  // Tag defintion:
  typedef Tag_true                                   Has_left_category;
  typedef Tag_true                                   Has_merge_category;
  typedef Tag_false                                  Has_do_intersect_category;

  typedef Arr_oblivious_side_tag                     Left_side_category;
  typedef Arr_oblivious_side_tag                     Bottom_side_category;
  typedef Arr_oblivious_side_tag                     Top_side_category;
  typedef Arr_oblivious_side_tag                     Right_side_category;

private:
  typedef Arr_polyline_traits_2<Segment_traits_2>    Self;

  // Data members:
  Segment_traits_2   m_seg_traits;           // The base segment-traits class.

private:
  enum { INVALID_INDEX = 0xffffffff };
  
public:
  /*! Default constructor */
  Arr_polyline_traits_2() : m_seg_traits() {}

  /*! Obtain the segment traits.
   * \return the segment traits.
   */
  const Segment_traits_2* segment_traits_2() const { return &m_seg_traits; }

  /// \name Types and functors inherited from the base segment traits.
  //@{

  // Traits types:
  typedef typename Segment_traits_2::Point_2            Point_2;
  typedef typename Segment_traits_2::Curve_2            Segment_2;

  typedef _Polyline_2<Segment_traits_2>                 Curve_2;
  typedef _X_monotone_polyline_2<Segment_traits_2>      X_monotone_curve_2;

  typedef typename Segment_traits_2::Multiplicity       Multiplicity;

  /*! Compare the x-coordinates of two points. */
  typedef typename Segment_traits_2::Compare_x_2        Compare_x_2;

  /*! Get a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object() const
  {
    return m_seg_traits.compare_x_2_object();
  }

  /*! Compare two points lexigoraphically: by x, then by y. */
  typedef typename Segment_traits_2::Compare_xy_2       Compare_xy_2;
  
  /*! Get a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const
  {
    return m_seg_traits.compare_xy_2_object();
  }

  ///@}

  /// \name Basic predicate functors(based on the segment traits).
  //@{

  class Construct_min_vertex_2 {
  private:
    const Segment_traits_2  * m_seg_traits;

  public:

    /*! Constructor. */
    Construct_min_vertex_2(const Segment_traits_2 * traits) : m_seg_traits(traits)
    {}

    /*!
     * Get the left endpoint of the x-monotone curve(segment).
     * \param cv The polyline curve.
     * \return The left endpoint.
     */
    const Point_2 operator()(const X_monotone_curve_2 & cv) const
    {
      CGAL_assertion(cv.size() > 0);

      return m_seg_traits->construct_min_vertex_2_object()(cv[0]);
    }
  };
    
  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  {
    return Construct_min_vertex_2(&m_seg_traits);
  }

  class Construct_max_vertex_2 {
  private:
    const Segment_traits_2  * m_seg_traits;

  public:
    /*! Constructor. */
    Construct_max_vertex_2(const Segment_traits_2 * traits) : m_seg_traits(traits)
    {}

    /*!
     * Get the right endpoint of the x-monotone curve(segment).
     * \param cv The polylinecurve.
     * \return The right endpoint.
     */
    const Point_2 operator()(const X_monotone_curve_2 & cv) const
    {
      CGAL_assertion(cv.size() > 0);

      return m_seg_traits->construct_max_vertex_2_object()(cv[cv.size() - 1]);
    }
  };

  /*! Get a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  {
    return Construct_max_vertex_2(&m_seg_traits);
  }

  class Is_vertical_2 {
  private:
    const Segment_traits_2  * m_seg_traits;

  public:
    /*! Constructor. */
    Is_vertical_2(const Segment_traits_2 * traits) : m_seg_traits(traits) {}

    /*!
     * Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2 & cv) const
    {
      // An x-monotone polyline can represent a vertical segment only if it
      // is comprised of vertical segments. If the first segment is vertical,
      // all segments are vertical in an x-monotone polyline
      return (m_seg_traits->is_vertical_2_object()(cv[0]));
    }
  };
  
  /*! Get an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const
  {
    return Is_vertical_2(&m_seg_traits);
  }

  class Compare_y_at_x_2 {
  private:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Self;
    const Segment_traits_2  * m_seg_traits;

  public:

    /*! Constructor. */
    Compare_y_at_x_2(const Segment_traits_2 * traits) : m_seg_traits(traits) {}

    /*!
     * Return the location of the given point with respect to the input curve.
     * \param cv The polyline curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator()(const Point_2 & p,
                                 const X_monotone_curve_2 & cv) const
    {
      // Get the index of the segment in cv containing p.
      unsigned int i = Self::_locate(m_seg_traits, cv, p);
      CGAL_precondition(i != INVALID_INDEX);

      // Compare the segment cv[i] and p.
      return m_seg_traits->compare_y_at_x_2_object()(p, cv[i]);
    }
  };
  friend class Compare_y_at_x_2;

  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  {
    return Compare_y_at_x_2(&m_seg_traits);
  }

  class Compare_y_at_x_left_2 {
  private:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Self;

    const Segment_traits_2  * m_seg_traits;

  public:
    /*! Constructor. */
    Compare_y_at_x_left_2(const Segment_traits_2 * traits) : m_seg_traits(traits)
    {}

    /*!
     * Compare the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param cv1 The first polyline curve.
     * \param cv2 The second polyline curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined(lexicographically) to its left.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2 & cv1,
                                 const X_monotone_curve_2 & cv2,
                                 const Point_2 & p) const
    {
      // Get the indices of the segments in cv1 and cv2 containing p and
      // defined to its left.
      unsigned int i1 = Self::_locate_side(m_seg_traits, cv1, p, false);
      unsigned int i2 = Self::_locate_side(m_seg_traits, cv2, p, false);

      CGAL_precondition(i1 != INVALID_INDEX);
      CGAL_precondition(i2 != INVALID_INDEX);

      // Compare cv1[i1] and cv2[i2] at p's left.
      return m_seg_traits->compare_y_at_x_left_2_object()(cv1[i1], cv2[i2], p);
    }
  };
  friend class Compare_y_at_x_left_2;

  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  {
    return Compare_y_at_x_left_2(&m_seg_traits);
  }

  class Compare_y_at_x_right_2 {
  private:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Self;
    const Segment_traits_2  * m_seg_traits;

  public:
    /*! Constructor. */
    Compare_y_at_x_right_2(const Segment_traits_2 * traits) :
      m_seg_traits(traits)
    {}

    /*!
     * Compare the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined(lexicographically) to its right.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2 & cv1,
                                 const X_monotone_curve_2 & cv2,
                                 const Point_2 & p) const
    {
      // Get the indices of the segments in cv1 and cv2 containing p and
      // defined to its right.
      unsigned int i1 = Self::_locate_side(m_seg_traits, cv1, p, true);
      unsigned int i2 = Self::_locate_side(m_seg_traits, cv2, p, true);

      CGAL_precondition(i1 != INVALID_INDEX);
      CGAL_precondition(i2 != INVALID_INDEX);

      // Compare cv1[i1] and cv2[i2] at p's right.
      return m_seg_traits->compare_y_at_x_right_2_object()(cv1[i1],cv2[i2], p);
    }
  };
  friend class Compare_y_at_x_right_2;

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  {
    return Compare_y_at_x_right_2(&m_seg_traits);
  }

  class Equal_2 {
  private: 
    const Self * m_poly_traits;

  public:

    /*! Constructor. */
    Equal_2(const Self * poly_tr) : m_poly_traits(poly_tr) {}

    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same;(false) otherwise.
     */
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      return m_poly_traits->m_seg_traits.equal_2_object()(p1, p2);
    }

    /*!
     * Check if the two x-monotone curves are the same(have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return(true) if the two curves are the same;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2 & cv1,
                    const X_monotone_curve_2 & cv2) const
    {
      unsigned int n1 = cv1.size();
      unsigned int n2 = cv2.size();

      // Check the pairwise equality of the contained segments.    
      typename Segment_traits_2::Equal_2 equal =
        m_poly_traits->m_seg_traits.equal_2_object();
      typename Segment_traits_2::Compare_x_2 compare_x =
        m_poly_traits->m_seg_traits.compare_x_2_object(); 
      typename Segment_traits_2::Compare_y_at_x_2 compare_y_at_x =
        m_poly_traits->m_seg_traits.compare_y_at_x_2_object();
      typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
        m_poly_traits->m_seg_traits.construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        m_poly_traits->m_seg_traits.construct_max_vertex_2_object();
      Is_vertical_2 is_vertical = m_poly_traits->is_vertical_2_object();
      Point_2 point1,point2;
      Comparison_result res_x;
      Comparison_result res_y_at_x;
      unsigned int i = 0, j = 0;
      
      // the first and last points of the segments should be equal.
      bool res = equal(min_vertex(cv1[0]),min_vertex(cv2[0]));
      if (!res)
        return false;
      res = equal(max_vertex(cv1[n1-1]),max_vertex(cv2[n2-1]));
      if (!res)
        return false;

      // If the first and last points are equal and the curves are vertical,
      // it means that it is equal.
      bool ver1 = is_vertical(cv1);
      bool ver2 = is_vertical(cv2);
      // both curves are vertical and therefore equal. 
      if (ver1 && ver2)
        return true;
      // one is vertical and the other is not - hence not equal.
      if ((ver1 && !ver2) || (ver2 && !ver1))
        return false;
      
      // If we arrived here it means that the first and last point of the
      // curve are equal.
      while ((i < n1-1) || (j < n2-1)) {
        point1 = max_vertex(cv1[i]);
        point2 = max_vertex(cv2[j]);

        res = equal(point1, point2);
        // Easy case - the two points are equal
        if (res) {
          ++i;
          ++j;
        }
        else {
          res_x = compare_x(point1,point2);
          // Check if the different point is a collinear point situated on
          // the line between its two neighbors.
          if (res_x == SMALLER) { 
            res_y_at_x = compare_y_at_x(point1,cv2[j]);
            if (res_y_at_x == EQUAL)
              ++i;
            else
              return false;
          }
          else if(res_x == LARGER) {
            res_y_at_x = compare_y_at_x(point2,cv1[i]);
            if (res_y_at_x == EQUAL)
              ++j;
            else
              return false;
          }
          else {
            return false;
          }
        }
      }
      return true;
    }
  };
 
  /*! Get an Equal_2 functor object. */
  
  Equal_2 equal_2_object() const
  {
    return Equal_2(this);
  }
  ///@}

  /// \name Construction functors(based on the segment traits).
  //@{

  class Make_x_monotone_2 {
  private:
    const Segment_traits_2  * m_seg_traits;

  public:
    /*! Constructor. */
    Make_x_monotone_2(const Segment_traits_2 * traits) : m_seg_traits(traits) {}

    /*!
     * Cut the given curve into x-monotone subcurves and insert them into the
     * given output iterator. As segments are always x_monotone, only one
     * object will be contained in the iterator.
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. The output
     *           object is a wrapper of either an X_monotone_curve_2, or - in
     *           case the input segment is degenerate - a Point_2 object.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator()(const Curve_2 & cv, OutputIterator oi) const
    { 
      // Go over all points in the input curve.
      typename Curve_2::const_iterator       ps = cv.begin();
      typename Curve_2::const_iterator       end = cv.end();

      // Empty polyline:
      if (ps == end)
        return oi;

      typename Curve_2::const_iterator       pt = ps;
      ++pt;

      if (pt == end) {
        // The polyline contains a single isolated point:
        *oi++ = make_object(*ps);
        return oi;
      }

      // Locate points where the x-order changes:
      typename Segment_traits_2::Compare_x_2     compare_x =
        m_seg_traits->compare_x_2_object();
      typename Segment_traits_2::Compare_xy_2    compare_xy =
        m_seg_traits->compare_xy_2_object();
      Comparison_result                          x_res;
      Comparison_result                          xy_res;
      typename Curve_2::const_iterator           x_begin = ps;
      Comparison_result                          curr_x_res;
      Comparison_result                          curr_xy_res;
     
      x_res = compare_x(*ps, *pt);
      if (x_res != EQUAL)
        xy_res = x_res;
      else
        xy_res = compare_xy(*ps, *pt);

      ++ps; ++pt;
      while (pt != end) {
        curr_x_res = compare_x(*ps, *pt);
        if (curr_x_res != EQUAL)
          curr_xy_res = curr_x_res;
        else
          curr_xy_res = compare_xy(*ps, *pt);

        if (curr_x_res != x_res || curr_xy_res != xy_res) {
          // Create a new x-monotone polyline from the range of points
          // [x_begin, pt):
          *oi++ = make_object(X_monotone_curve_2(x_begin, pt));

          x_begin = ps;
          x_res = curr_x_res;
          xy_res = curr_xy_res;
        }

        ++ps; ++pt;
      }

      // Create an x-monotone polyline from the remaining points.
      CGAL_assertion(x_begin != end);
      *oi++ = make_object(X_monotone_curve_2(x_begin, end));
      return oi;
    }
  };
  
  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  {
    return Make_x_monotone_2(&m_seg_traits);
  }

  class Split_2 
  {
  private:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Self;

    const Segment_traits_2 * m_seg_traits;

  public:
    /*! Constructor. */
    Split_2(const Segment_traits_2 * traits) : m_seg_traits(traits) {}

  public:
    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve(p is its right endpoint).
     * \param c2 Output: The right resulting subcurve(p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2 & cv, const Point_2 & p,
                    X_monotone_curve_2 & c1, X_monotone_curve_2 & c2) const
    {
      typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
        m_seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        m_seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Equal_2 equal = m_seg_traits->equal_2_object();

      // Make sure the split point is not one of the curve endpoints.
      CGAL_precondition(!equal(min_vertex(cv[0]), p));
      CGAL_precondition(!equal(max_vertex(cv[cv.size() - 1]), p));

      // Locate the segment on the polyline cv that contains p.
      unsigned int i = Self::_locate(m_seg_traits, cv, p);
      CGAL_precondition(i != INVALID_INDEX);

      // Clear the output curves.
      c1.clear(); 
      c2.clear();

      // Push all segments labeled(0, 1, ... , i-1) into c1.
      unsigned int j;
      for (j = 0; j < i; ++j)
        c1.push_back(cv[j]);

      // Check whether the split point is cv[i]'s source of target.
      if (equal(max_vertex(cv[i]), p)) {
        // The entire i'th segment belongs to c1:
        c1.push_back(cv[i]);
      } else if (equal(min_vertex(cv[i]), p)) {
        // The entire i'th segments belongs to c2:
        c2.push_back(cv[i]);
      } else {
        // The i'th segment should be split: The left part(seg1) goes to cv1,
        // and the right part(seg2) goes to cv2.
        Segment_2   seg1, seg2;
        m_seg_traits->split_2_object()(cv[i], p, seg1, seg2);

        c1.push_back(seg1);
        c2.push_back(seg2);
      }

      // Push all segments labeled(i+1, i+2, ... , n-1) into cv1.
      unsigned int n = cv.size();

      for (j = i+1; j < n; ++j)
        c2.push_back(cv[j]);
    }
  };
  friend class Split_2;

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object() const
  {
    return Split_2(&m_seg_traits);
  }

  class Intersect_2 {
  private:
    typedef Arr_polyline_traits_2<Segment_traits_2>       Self;
    const Segment_traits_2 * m_seg_traits;

  public:
    /*! Constructor. */
    Intersect_2(const Segment_traits_2 * traits) : m_seg_traits(traits) {}

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
    OutputIterator operator()(const X_monotone_curve_2 & cv1,
                              const X_monotone_curve_2 & cv2,
                              OutputIterator oi)
    {
      typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
        m_seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        m_seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Equal_2                equal =
        m_seg_traits->equal_2_object();
      typename Segment_traits_2::Compare_xy_2           compare_xy =
        m_seg_traits->compare_xy_2_object();
      typename Segment_traits_2::Intersect_2            intersect =
        m_seg_traits->intersect_2_object();
      typename Segment_traits_2::Compare_y_at_x_2       compare_y_at_x =
        m_seg_traits->compare_y_at_x_2_object();
      
      const unsigned int n1 = cv1.size();
      const unsigned int n2 = cv2.size();
      unsigned int       i1 = 0;
      unsigned int       i2 = 0;
      X_monotone_curve_2 ocv;           // Used to represent overlaps.

      Comparison_result left_res = compare_xy(min_vertex(cv1[i1]),
                                              min_vertex(cv2[i2]));
      
      if (left_res == SMALLER) {
        // cv1's left endpoint is to the left of cv2's left endpoint:
        // Locate the index i1 of the segment in cv1 which contains cv2's
        // left endpoint.
        i1 = Self::_locate(m_seg_traits, cv1, min_vertex(cv2[i2]));
        if (i1 == INVALID_INDEX)
          return oi;

        if (equal(max_vertex(cv1[i1]), min_vertex(cv2[i2]))) {
          ++i1;
          left_res = EQUAL;
        }
      }  else if (left_res == LARGER) {
        // cv1's left endpoint is to the right of cv2's left endpoint:
        // Locate the index i2 of the segment in cv2 which contains cv1's
        // left endpoint.
        i2 = Self::_locate(m_seg_traits, cv2, min_vertex(cv1[i1]));

        if (i2 == INVALID_INDEX)
          return oi;

        if (equal(max_vertex(cv2[i2]), min_vertex(cv1[i1]))) {
          ++i2;
          left_res = EQUAL;
        }
      }

      // Check if the the left endpoint lies on the other polyline.
      bool left_coincides = (left_res == EQUAL);
      bool left_overlap = false;

      if (left_res == SMALLER)
        left_coincides =(compare_y_at_x(min_vertex(cv2[i2]), cv1[i1]) == EQUAL);
      else if (left_res == LARGER)
        left_coincides =(compare_y_at_x(min_vertex(cv1[i1]), cv2[i2]) == EQUAL);

      // The main loop: Go simultaneously over both polylines.
      Comparison_result right_res = left_res;
      bool right_coincides = left_coincides;
      bool right_overlap = false;

      while (i1 < n1 && i2 < n2) {
        right_res = compare_xy(max_vertex(cv1[i1]), max_vertex(cv2[i2]));

        right_coincides = (right_res == EQUAL);
        if (right_res == SMALLER)
          right_coincides = (compare_y_at_x(max_vertex(cv1[i1]),
                                            cv2[i2]) == EQUAL);
        else if (right_res == LARGER)
          right_coincides = (compare_y_at_x(max_vertex(cv2[i2]),
                                            cv1[i1]) == EQUAL);

        right_overlap = false;
        
        if (!right_coincides && !left_coincides) {
          // Non of the endpoints of the current segment of one polyline
          // coincides with the curent segment of the other polyline:
          // Output the intersection if exists.
          oi = intersect(cv1[i1], cv2[i2], oi);
        } else if (right_coincides && left_coincides) {
          // An overlap exists between the current segments of the polylines:
          // Output the overlapping segment.
          right_overlap = true;
          if (left_res == SMALLER) {
            if (right_res == SMALLER) {
              Segment_2 seg(min_vertex(cv2[i2]), max_vertex(cv1[i1]));
              ocv.push_back(seg);
            } else {
              Segment_2 seg(min_vertex(cv2[i2]), max_vertex(cv2[i2]));
              ocv.push_back(seg);
            }
          } else {
            if (right_res == SMALLER) {
              Segment_2 seg(min_vertex(cv1[i1]), max_vertex(cv1[i1]));
              ocv.push_back(seg);
            } else {
              Segment_2 seg(min_vertex(cv1[i1]), max_vertex(cv2[i2]));
              ocv.push_back(seg);
            }
          }
        } else if (left_coincides && !right_coincides) {
          // The left point of the current segment of one polyline
          // coincides with the current segment of the other polyline.
          if (left_overlap) {
            // An overlap occured at the previous iteration:
            // Output the overlapping polyline.
            CGAL_assertion(ocv.size() > 0);
            *oi++ = make_object(ocv);
            ocv.clear();
          } else {
            // The left point of the current segment of one polyline
            // coincides with the current segment of the other polyline, and
            // no overlap occured at the previous iteration:
            // Output the intersection point. The derivative of at least one of
            // the polylines is not defined at this point, so we give it
            // multiplicity 0.
            if (left_res == SMALLER) {
              std::pair<Point_2,Multiplicity>  p(min_vertex(cv2[i2]), 0);
              *oi++ = make_object(p);
            } else {
              std::pair<Point_2,Multiplicity>  p(min_vertex(cv1[i1]), 0);
              *oi++ = make_object(p);
            }
          }
        }

        // Proceed forward.
        if (right_res != SMALLER)
          ++i2;
        if (right_res != LARGER)
          ++i1;

        left_res = (right_res == SMALLER) ? LARGER :
          (right_res == LARGER) ? SMALLER : EQUAL;
        
        left_coincides = right_coincides;
        left_overlap = right_overlap;
      }

#if 0
      std::cout << "left coincides: " << left_coincides << std::endl;
      std::cout << "right coincides: " << right_coincides << std::endl;
      std::cout << "right res: " << right_res << std::endl;
      std::cout << "left res: " << left_res << std::endl;
#endif
      
      // Output the remaining overlapping polyline, if necessary.
      if (ocv.size() > 0) {
        *oi++ = make_object(ocv);
      } else if (right_coincides) {
        if (right_res == SMALLER) {
          std::pair<Point_2,Multiplicity> ip(max_vertex(cv1[i1 - 1]), 0);
          *oi++ = make_object(ip);
        } else if (right_res == LARGER) {
          std::pair<Point_2,Multiplicity> ip(max_vertex(cv2[i2 - 1]), 0);
          *oi++ = make_object(ip);
        } else if (i1 > 0) {
          std::pair<Point_2,Multiplicity> ip(max_vertex(cv1[i1 - 1]), 0);
          *oi++ = make_object(ip);
        } else {
          CGAL_assertion(i2 > 0);
          std::pair<Point_2,Multiplicity> ip(max_vertex(cv2[i2 - 1]), 0);
          *oi++ = make_object(ip);
        }
      }
      return oi;
    }
  };
  friend class Intersect_2;
  
  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const
  {
    return Intersect_2(&m_seg_traits);
  }

  class Are_mergeable_2 {
  private:
    const Segment_traits_2 * m_seg_traits;

  public:
    /*! Constructor. */
    Are_mergeable_2(const Segment_traits_2 * traits) : m_seg_traits(traits) {}

    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return(true) if the two curves are mergeable, that is, they share a
     * common endpoint;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2 & cv1,
                    const X_monotone_curve_2 & cv2) const
    {
      typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
        m_seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        m_seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Equal_2 equal = m_seg_traits->equal_2_object();

      typename Segment_traits_2::Is_vertical_2 is_vertical =
        m_seg_traits->is_vertical_2_object();
      
      const unsigned int n1 = cv1.size();
      const unsigned int n2 = cv2.size();

      bool ver1 = is_vertical(cv1[0]);
      bool ver2 = is_vertical(cv2[0]);

      return ((equal(max_vertex(cv1[n1 - 1]), min_vertex(cv2[0])) ||
               equal(max_vertex(cv2[n2 - 1]), min_vertex(cv1[0]))) &&
              ((ver1 && ver2) || (!ver1 && !ver2)));
    }
  };
  
  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object() const
  {
    return Are_mergeable_2(&m_seg_traits);
  }

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>     Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;
    
    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_polyline_traits_2<Segment_traits_2>;
    
  public:
    /*!
     * Merge two given x-monotone curves into a single curve(segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable.
     */
    void operator()(const X_monotone_curve_2 & cv1,
                    const X_monotone_curve_2 & cv2,
                    X_monotone_curve_2 & c) const
    {
      CGAL_precondition(m_traits->are_mergeable_2_object()(cv1, cv2));

      const Segment_traits_2* seg_traits = m_traits->segment_traits_2();

      Construct_min_vertex_2 min_vertex =
        m_traits->construct_min_vertex_2_object();
      Construct_max_vertex_2 max_vertex =
        m_traits->construct_max_vertex_2_object();
      Equal_2 equal = m_traits->equal_2_object();
      
      const unsigned int n1 = cv1.size();
      const unsigned int n2 = cv2.size();
      unsigned int       i;

      c.clear();
      if (equal(max_vertex(cv1), min_vertex(cv2))) {
        // cv2 extends cv1 to the right:
        for (i = 0; i < n1 - 1; ++i)
          c.push_back(cv1[i]);

        // Try to merge tthe to contiguous line segments:
        if (seg_traits->are_mergeable_2_object()(cv1[n1 - 1], cv2[0])) {
          Segment_2       seg;
          seg_traits->merge_2_object()(cv1[n1 - 1], cv2[0], seg);
          c.push_back(seg);
        } else {
          c.push_back(cv1[n1 - 1]);
          c.push_back(cv2[0]);
        }

        for (i = 1; i < n2; ++i)
          c.push_back(cv2[i]);
      } else {
        CGAL_precondition(equal(max_vertex(cv2), min_vertex(cv1)));
        
        // cv1 extends cv2 to the right:
        for (i = 0; i < n2 - 1; ++i)
          c.push_back(cv2[i]);

        // Try to merge tthe to contiguous line segments:
        if (seg_traits->are_mergeable_2_object()(cv2[n2 - 1], cv1[0])) {
          Segment_2       seg;
          seg_traits->merge_2_object()(cv2[n2 - 1], cv1[0], seg);
          c.push_back(seg);
        } else {
          c.push_back(cv2[n2 - 1]);
          c.push_back(cv1[0]);
        }

        for (i = 1; i < n1; ++i)
          c.push_back(cv1[i]);
      }
    }
  };
  
  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object() const { return Merge_2(this); }
  ///@}
  
  /// \name Functor definitions for the landmarks point-location strategy.
  //@{
  typedef typename Segment_traits_2::Approximate_number_type  
  Approximate_number_type;
  typedef typename Segment_traits_2::Approximate_2    Approximate_2;


  /*! Get an Approximate_2 functor object. */
  Approximate_2 approximate_2_object () const
  {
    return m_seg_traits.approximate_2_object();
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
      // Construct a polyline containing just two points:
      Point_2   pts[2];

      pts[0] = p; pts[1] = q;
      return (X_monotone_curve_2 (pts + 0, pts + 2));
    }
  };

  /*! Get a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object () const
  {
    return Construct_x_monotone_curve_2();
  }
  //@}

private:
  /*!
   * Return the index of the segment in the polyline that contains the
   * point q in its x-range. The function performs a binary search, so if the
   * point q is in the x-range of the polyline with n segments, the segment
   * containing it can be located in O(log n) operations.
   * \param cv The polyline curve.
   * \param q The point.
   * \return An index i such that q is in the x-range of cv[i].
   *         If q is not in the x-range of cv, returns INVALID_INDEX.
   */
  static unsigned int _locate(const Segment_traits_2 * m_seg_traits,
                              const X_monotone_curve_2 & cv,
                              const Point_2 & q)
  {
    typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
      m_seg_traits->construct_min_vertex_2_object();
    typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
      m_seg_traits->construct_max_vertex_2_object();

    unsigned int from = 0;
    unsigned int to = cv.size() - 1;

    if (m_seg_traits->is_vertical_2_object()(cv[0])) {
      typename Segment_traits_2::Compare_xy_2 compare_xy =
        m_seg_traits->compare_xy_2_object();

      // First check whether the polyline curve really contains q in its
      // xy-range:

      Comparison_result res_from = compare_xy(min_vertex(cv[from]), q);    
      if (res_from == EQUAL) return from;
    
      Comparison_result res_to = compare_xy(max_vertex(cv[to]), q);
      if (res_to == EQUAL) return to;
    
      typename Segment_traits_2::Compare_x_2 compare_x =
        m_seg_traits->compare_x_2_object();

      Comparison_result res = compare_x(max_vertex(cv[to]), q);
      if (res != EQUAL) return INVALID_INDEX;

      //// q is not in the x-range of cv:
      //if (res_from == res_to) return INVALID_INDEX;
      
      // Perform a binary search to locate the segment that contains q in its
      // xy-range:
      while (to > from) {
        unsigned int mid = (from + to) / 2;
        if (mid > from) {
          Comparison_result res_mid = compare_xy(min_vertex(cv[mid]), q);
          if (res_mid == EQUAL) return mid;
          if (res_mid == res_from) from = mid;
          else to = mid - 1;
        } else {
          CGAL_assertion(mid < to);
          Comparison_result res_mid = compare_xy(max_vertex(cv[mid]), q);
          if (res_mid == EQUAL) return mid;
          if (res_mid == res_to) to = mid;
          else from = mid + 1;
        }
      }
      // In case (from == to), and we know that the polyline contains the q:
      CGAL_assertion(from == to);
      return from;
    }
    
    typename Segment_traits_2::Compare_x_2 compare_x =
      m_seg_traits->compare_x_2_object();

    // First check whether the polyline curve really contains q in its x-range.
    Comparison_result res_from = compare_x(min_vertex(cv[from]), q);    
    if (res_from == EQUAL) return from;
    
    Comparison_result res_to = compare_x(max_vertex(cv[to]), q);
    if (res_to == EQUAL) return to;
    
    // q is not in the x-range of cv:
    if (res_from == res_to) return INVALID_INDEX;

    // Perform a binary search and locate the segment that contains q in its
    // x-range:
    while (to > from) {
      unsigned int mid = (from + to) / 2;

      if (mid > from) {
        Comparison_result res_mid = compare_x(min_vertex(cv[mid]), q);
        if (res_mid == EQUAL) return mid;
        if (res_mid == res_from) from = mid;
        else to = mid - 1;
      } else {
        CGAL_assertion(mid < to);
        Comparison_result res_mid = compare_x(max_vertex(cv[mid]), q);
        if (res_mid == EQUAL) return mid;
        if (res_mid == res_to) to = mid;
        else from = mid + 1;
      }
    }

    // In case(from == to), and we know that the polyline contains the q:
    CGAL_assertion(from == to);
    return from;
  }

  /*!
   * Find the index of the segment in the polyline that is defined to the
   * left(or to the right) of the point q.
   * \param cv The polyline curve.
   * \param q The point.
   * \param to_right(true) if we wish to locate a segment to the right of q,
   *               (false) if we wish to locate a segment to its right.
   * \return An index i such that segments[i] is defined to the left(or to the
   *         right) of q, or INVALID_INDEX if no such segment exists.
   */
  static unsigned int _locate_side(const Segment_traits_2 * m_seg_traits,
                                   const X_monotone_curve_2 & cv,
                                   const Point_2 & q, const bool & to_right)
  {
    // First locate a segment segments[i] that contains q in its x-range.
    unsigned int i = _locate(m_seg_traits, cv, q);
    if (i == INVALID_INDEX)
      return INVALID_INDEX;
   
    typename Segment_traits_2::Equal_2 equal = m_seg_traits->equal_2_object();

    if (equal(m_seg_traits->construct_min_vertex_2_object()(cv[i]), q)) {
      // q is the left endpoint of the i'th segment:
      if (to_right)
        return i;
      else if (i == 0)
        return INVALID_INDEX;
      else
        return i - 1;
    }

    if (equal(m_seg_traits->construct_max_vertex_2_object()(cv[i]), q)) {
      // q is the right endpoint of the i'th segment:
      if (!to_right)
        return i;
      else if (i == (cv.size() - 1))
        return INVALID_INDEX;
      else
        return i + 1;
    }
    
    // In case q is in cv[i]'s interior:
    return i;
  }
};


} //namespace CGAL

#endif
