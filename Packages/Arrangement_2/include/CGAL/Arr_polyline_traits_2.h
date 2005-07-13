// Copyright (c) 2003  Tel-Aviv University (Israel).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POLYLINE_TRAITS_2_H
#define CGAL_ARR_POLYLINE_TRAITS_2_H

/*! \file
 * The traits-class for the linear piece-wiese (polyline) type of curves of the
 * arrangement package.
 */

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/intersections.h>

#include <list>
#include <fstream>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class T_Segment_traits> class Polyline_2;

template <class T_Segment_traits>
class Arr_polyline_traits_2 {
public:
  typedef Tag_true                                      Has_left_category;
  typedef Tag_true                                      Has_reflect_category;  
  typedef Tag_true                                      Has_merge_category;
    
  typedef T_Segment_traits                              Segment_traits_2;
  typedef typename Segment_traits_2::Kernel             Kernel;

protected:
  typedef Arr_polyline_traits_2<Segment_traits_2>       Self;

public:
  /*! Default constructor */
  Arr_polyline_traits_2() {}

  /// \name Types and functor inherited from the segment traits
  //@{

  // Traits types:
  typedef typename Segment_traits_2::Point_2            Point_2;
  typedef typename Segment_traits_2::Curve_2            Segment_2;

  typedef Polyline_2<Segment_traits_2>                  Curve_2;
  typedef Polyline_2<Segment_traits_2>                  X_monotone_curve_2;

  /*! Compare the x-coordinates of two points */
  typedef typename Segment_traits_2::Compare_x_2        Compare_x_2;

  /*! Get a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object() const
  {
    return Compare_x_2();
  }

  /*! Compare two points lexigoraphically; by x, then by y */
  typedef typename Segment_traits_2::Compare_xy_2       Compare_xy_2;
  
  /*! Get a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const
  {
    return Compare_xy_2();
  }

  ///@}

  /// \name Basic functors introduced here (based on the segment traits)
  //@{

  class Construct_min_vertex_2 {
  public:
    /*! Get the left endpoint of the x-monotone curve(segment).
     * \param cv The polyline curve.
     * \return The left endpoint.
     */
    const Point_2 & operator()(const X_monotone_curve_2 & cv) const
    {
      CGAL_assertion(cv.size() > 0);
      Segment_traits_2 seg_traits;
      return seg_traits.construct_min_vertex_2_object()(cv[0]);
    }
  };
    
  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  {
    return Construct_min_vertex_2();
  }

  class Construct_max_vertex_2 {
  public:
    /*! Get the right endpoint of the x-monotone curve(segment).
     * \param cv The polylinecurve.
     * \return The right endpoint.
     */
    const Point_2 & operator()(const X_monotone_curve_2 & cv) const
    {
      CGAL_assertion(cv.size() > 0);
      Segment_traits_2 seg_traits;
      return seg_traits.construct_max_vertex_2_object()(cv[cv.size() - 1]);
    }
  };

  /*! Get a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  {
    return Construct_max_vertex_2();
  }

  class Is_vertical_2 {
  public:
    /*! Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2 & cv) const
    {
      // An x-monotone polyline can represent a vertical segment only if it
      // is comprised of a single vertical segment.
      Segment_traits_2 seg_traits;
      return (cv.size() == 1 && seg_traits.is_vertical_object()(cv[0]));
    }
  };
  
  /*! Get an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const
  {
    return Is_vertical_2();
  }

  class Compare_y_at_x_2 {
  public:
    /*! Return the location of the given point with respect to the input curve.
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
      unsigned int i = locate_point(cv, p);
      CGAL_precondition(i != 0xffffffff);

      // Compare the segment cv[i] and p.
      Segment_traits_2 seg_traits;
      return seg_traits.compare_y_at_x_2_object()(p, cv[i]);
    }
  };
  
  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  {
    return Compare_y_at_x_2();
  }

  class Equal_2 {
  public:
    /*! Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same;(false) otherwise.
     */
    bool operator()(const Point_2 & p1, const Point_2 & p2) const
    {
      Segment_traits_2 seg_traits;
      return seg_traits.equal_2_object()(p1, p2);
    }

    /*! Check if the two x-monotone curves are the same(have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2 & cv1,
                    const X_monotone_curve_2 & cv2) const
    {
      // The two curves must contain the same number of segments.
      unsigned int n1 = cv1.size();
      unsigned int n2 = cv2.size();
      if (n1 != n2) return false;

      // Check the pairwise equality of the contained segments. 
      Segment_traits_2 seg_traits;
      Equal_2 equal = seg_traits.equal_2_object();
      for (unsigned int i = 0; i < n1; ++i)
        if (!equal(cv1[i], cv2[i])) return false;
      return true;
    }
  };

  /*! Get an Equal_2 functor object. */
  Equal_2 equal_2_object() const
  {
    return Equal_2();
  }

  class Compare_y_at_x_left_2 {
  public:
    /*! Compare the y value of two x-monotone curves immediately to the left
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
      // Get the indices of the segments in cv1 and cv2 containing q and
      // defined to its left.
      unsigned int i1 = locate_point_side(cv1, q, false);
      unsigned int i2 = locate_point_side(cv2, q, false);

      CGAL_precondition(i1 != 0xffffffff);
      CGAL_precondition(i2 != 0xffffffff);

      // Compare cv1[i1] and cv2[i2] at q's left.
      Segment_traits_2 seg_traits;
      return seg_traits.compare_y_at_x_left_object()(cv1[i1], cv2[i2], q);
    }
  };
  
  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  {
    return Compare_y_at_x_left_2();
  }

  class Compare_y_at_x_right_2 {
  public:
    /*! Compare the y value of two x-monotone curves immediately to the right
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
      // Get the indices of the segments in cv1 and cv2 containing q and
      // defined to its right.
      unsigned int i1 = locate_point_side(cv1, q, true);
      unsigned int i2 = locate_point_side(cv2, q, true);

      CGAL_precondition(i1 != 0xffffffff);
      CGAL_precondition(i2 != 0xffffffff);

      // Compare cv1[i1] and cv2[i2] at q's right.
      Segment_traits_2 seg_traits;
      return seg_traits.compare_y_at_x_right_2_object()(cv1[i1], cv2[i2], q);
    }
  };
  
  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  {
    return Compare_y_at_x_right_2();
  }

  class Reflect_xy_2 {
  public:
    /*! Reflect the given point through the origin.
     * \param pt The point to be reflected.
     * \return The reflected point.
     */
    Point_2 operator()(const Point_2 & pt) const
    {
      Segment_traits_2 seg_traits;
      return seg_traits.reflect_xy_2_object()(pt);
    }
    
    /*! Reflect the given segment through the origin.
     * \param cv The segment to be reflected.
     * \return The reflected segment.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2 & cv) const
    {
      Segment_traits_2 seg_traits;
      Reflect_xy_2 reflect = seg_traits.reflect_xy_2_object();
      for (unsigned int i = 0; i < cv.size(); ++i) reflect(cv[i]);
      cv.reverse();
    }
  };

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Reflect_xy_2 reflect_xy_2_object() const
  {
    return Reflect_xy_2();
  }
  
  ///@}

  /// \name functors introduced here(based on the segment traits)
  //@{

  class Make_x_monotone_2 {
  public:
    /*! Cut the given curve into x-monotone subcurves and insert them into the
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
      Segment_traits_2 seg_traits;
      Construct_min_vertex_2 min_ver =
        seg_traits.construct_min_vertex_2_object();
      Construct_max_vertex_2 max_ver =
        seg_traits.construct_max_vertex_2_object();
      Is_vertical_2 is_vertical = seg_traits.is_vertical_2_object();
      Equal_2 equal = seg_traits.equal_2_object();

      unsigned int n = curve.size();
      X_monotone_curve_2 x_cv;

      // Go over all curve segments.
      for (unsigned int i = 0; i < n; ++i) {
        if (is_vertical(curve[i])) {
          // The segment is vertical:
          if (x_cv.size() != 0) {
            // Cut the previous chain.
            *o++ = x_cv;
            x_cv.clear();
          }
          // Insert the vertical segment as a singleton polyline.
          x_cv.push_back(curve[i]);
          *o++ = x_cv;
          x_cv.clear();
          continue;
        }

        if (x_cv.size() == 0) {
          // There is no previous curve:
          // Just append the current curve to the current chain.
          x_cv.push_back(curve[i]);
          continue;
        }

        const Segment_traits_2::Point_2 & prev_max_vertex = max_ver(curve[i]-1);
        const Segment_traits_2::Point_2 & min_vertex = min_ver(curve[i]);
        if (!equal(prev_min_vertex, min_vertex)) {
          // The previous max is not identical to the current min:
          // Cut the previous chain.
          *o++ = x_cv;
          x_cv.clear();
          x_cv.push_back(curve[i]);
          continue;
        }
      }

      return o;
    }
  };
  
  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  {
    return Make_x_monotone_2();
  }

  class Split_2 {
  public:
    /*! Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve(p is its right endpoint).
     * \param c2 Output: The right resulting subcurve(p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2 & cv, const Point_2 & p,
                    X_monotone_curve_2 & c1, X_monotone_curve_2 & c2) const
    {
      // Locate the segment on the polyline cv that contains p.
      unsigned int i = locate_point(cv, p);
      CGAL_precondition(i != 0xffffffff);

      // Check preconditions.
      Segment_traits_2 seg_traits;
      Construct_min_vertex_2 min_ver =
        seg_traits.construct_min_vertex_2_object();
      Construct_max_vertex_2 max_ver =
        seg_traits.construct_max_vertex_2_object();
      Equal_2 equal = seg_traits.equal_2_object();
      CGAL_precondition(!equal(min_ver(cv), p));
      CGAL_precondition(!equal(max_ver(cv), p));

      // Clear the output curves.
      c1.clear(); 
      c2.clear();

      // Push all segments labeled(0, 1, ... , i-1) into c1.
      unsigned int j;
      for (j = 0; j < i; ++j) c1.push_back(cv[j]);

      // Check whether the split point is cv[i]'s source of target.
      if (equal(min_ver(cv[i]), p)) c2.push_back(cv[i]);
      else if (equal(max_ver(cv[i]), p)) c1.push_back(cv[i]);
      else {
        // cv[i] should be split.
        Segment_2 cvi_1, cvi_2;
        seg_traits.split_2_object()(cv[i], cvi_1, cvi_2, p);

        // The first part should go into c1 and the second into c2.
        c1.push_back(cvi_1);
        c2.push_back(cvi_2);
      }

      // Push all segments labeled(i+1, i+2, ... , n-1) into cv1.
      unsigned int n = cv.size();
      for (j = i+1; j < n; ++j) c2.push_back(cv[j]);
      return;
    }
  };
  
  /*! Get a Split_2 functor object. */
  Split_2 split_2_object() const
  {
    return Split_2();
  }

  class Intersect_2 {
  public:
    /*! Find the intersections of the two given curves and insert them into the
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
                              OutputIterator oi) const
    {
    }
  };
  
  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const
  {
    return Intersect_2();
  }

  class Are_mergeable_2 {
  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable - if they are supported
     *         by the same line and share a common endpoint;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2 & cv1,
                    const X_monotone_curve_2 & cv2) const
    {
    }
  };
  
  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object() const
  {
    return Are_mergeable_2();
  }

  class Merge_2 {
  public:
    /*! Merge two given x-monotone curves into a single curve(segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      same line and share a common endpoint.
     */
    void operator()(const X_monotone_curve_2 & cv1,
                    const X_monotone_curve_2 & cv2,
                    X_monotone_curve_2 & c) const
    {
    }
  };
  
  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object() const
  {
    return Merge_2();
  }

  ///@}
  
  /// \name functors required for the landmark point location

  typedef double                          Approximate_number_type;

  class Approximate_2
  {
  public:

    /*! Return an approximation of a point coordinate.
     * \param p The exact point.
     * \param i The coordinate index(either 0 or 1).
     * \pre i is either 0 or 1.
     * \return An approximation of p's x-coordinate(if i == 0), or an 
     *         approximation of p's y-coordinate(if i == 1).
     */
    Approximate_number_type operator()(const Point_2 & p, unsigned int i) const
    {
    }
  };
  
  /*! Get an Approximate_2 functor object. */
  Approximate_2 approximate_2_object() const
  {
    return Approximate_2();
  }

  ///@}

  /*! Check whether the given point is in the x-range of the given curve.
   * In out case, the curve is a segment [s, t], check whether x(s)<=x(q)<=x(t)
   * or whether x(t)<=x(q)<=x(s).
   * \param cv The curve.
   * \param q The point.
   * \return (true) if q is in the x-range of cv.
   */
  bool point_in_x_range(const X_monotone_curve_2& cv, const Point_2& q) const
  {
    return (locate_point(cv, q) >= 0);
  }

  /*! Flip a given curve.
   * \param cv The input curve.
   * \return The flipped curve.
   */
  X_monotone_curve_2 curve_opposite(const X_monotone_curve_2 & cv) const
  {
    // Copy the segments in revered order while flipping each one individually.
    int                n = cv.size();
    X_monotone_curve_2 flip_cv;
    for (unsigned int i = 0; i < n; ++i)
      flip_cv.push_back(seg_traits.curve_opposite(cv[n - i - 1]));
    return flip_cv; 
  }

  /*!
   * Find the nearest intersection of the two given curves to the right of 
   * a given reference point.
   * Nearest is defined as the lexicographically nearest point, not including 
   * the point reference point itself.
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
   * there is an overlapping subcurve, that contains the reference point in
   * its x-range, the function should return an X_monotone_curve_2 whose 
   * interior is strictly to the right of the reference point(that is, whose
   * left endpoint is the projection of the reference point onto the 
   * overlapping subcurve).
   * NOTE: When there is an overlap we will always return a SEGMENT(i.e.,
   *       a polyline with 2 points) even if the overlap is actually a longer
   *       polyline, as this is still sufficient for the arrangement.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \return An empty object if there is no intersection to the right of p.
   *         An object wrapping a Point_2 in case of a simple intersection.
   *         An object wrapping an X_monotone_curve_2 in case of an overlap.
   */
  Object nearest_intersection_to_right(const X_monotone_curve_2 & cv1,
                                       const X_monotone_curve_2 & cv2,
                                       const Point_2 & p) const
  {
    Point_2 ps [2];
    int res;

    res = _nearest_intersection_to_side(cv1, cv2, p,
                                        true,     // To the right of p.
                                        ps[0], ps[1]);

    // Check if there is no intersection.
    if (res == 0)
      return Object();
    
    // Check if the intersection is a single point.
    if (res == 1)
      return CGAL::make_object(ps[0]);

    // Return a polyline(which is actaully a single segment).
    return CGAL::make_object(X_monotone_curve_2(ps + 0, ps + 2));
  }

  /*!
   * Find the nearest intersection of the two given curves to the left of 
   * a given reference point.
   * Nearest is defined as the lexicographically nearest point, not including 
   * the point reference point itself.
   * If the intersection of the two curves is an X_monotone_curve_2, that is,
   * there is an overlapping subcurve, that contains the reference point in
   * its x-range, the function should return an X_monotone_curve_2 whose 
   * interior is strictly to the left of the reference point(that is, whose
   * right endpoint is the projection of the reference point onto the 
   * overlapping subcurve).
   * NOTE: When there is an overlap we will always return a SEGMENT(i.e.,
   *       a polyline with 2 points) even if the overlap is actually a longer
   *       polyline, as this is still sufficient for the arrangement.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \return An empty object if there is no intersection to the left of p.
   *         An object wrapping a Point_2 in case of a simple intersection.
   *         An object wrapping an X_monotone_curve_2 in case of an overlap.
   */
  Object nearest_intersection_to_left(const X_monotone_curve_2 & cv1,
                                      const X_monotone_curve_2 & cv2,
                                      const Point_2 & p) const
  {
    Point_2   ps [2];
    int       res;

    res = _nearest_intersection_to_side(cv1, cv2,
					 p,
					 false,     // To the left of p.
					 ps[0], ps[1]);

    // Check if there is no intersection.
    if (res == 0)
      return Object();
    
    // Check if the intersection is a single point.
    if (res == 1)
      return CGAL::make_object(ps[0]);

    // Return a polyline(which is actaully a single segment).
    return CGAL::make_object(X_monotone_curve_2(ps + 0, ps + 2));
  }

  /*!
   * Check whether the two given curves overlap.
   * \patam cv1 The first curve.
   * \patam cv2 The second curve.
   * \return (true) if the two curves overlap in a one-dimensional subcurve
   *(i.e., not in a finite number of points). Otherwise, if they have a finite
   * number of intersection points, or none at all, return (false).
   */
  bool curves_overlap(const X_monotone_curve_2 & cv1,
                      const X_monotone_curve_2 & cv2) const
  {
    // Get the leftmost point of cv1 and of cv2.
    const Point_2& left1 = is_curve_to_right(cv1) ? curve_source(cv1) : 
                                                     curve_target(cv1);
    const Point_2& left2 = is_curve_to_right(cv2) ? curve_source(cv2) : 
                                                     curve_target(cv2);

    // Pick the righmost point of the two to start with.
    Point_2 p = (seg_traits.compare_x(left1, left2) == LARGER) ? left1 : left2;

    // Try to find an intersection to the right of p.
    Point_2 p1, p2;
    int res;

    while ((res = _nearest_intersection_to_side(cv1, cv2, p,
                                                true,     // To the right of p.
                                                p1, p2)) != 0)
    {
      // Check if an overlap has been detected:
      if (res == 2)
	return true;

      // Otherwise, p1==p2 is the next intersection point to the right.
      // Continue checking the intersections to its right.
      p = p1;
    }

    // If we reached here, no overlap has been detected.
    return false;
  }

private:

  /*!
   * Return the index of the segment in the polyline cv that contains the
   * point q in its x-range. The function performs a binary search, so if the
   * point q is in the x-range of the polyline with n segments, the segment
   * containing it can be located in O(log n) operations.
   * \param cv The polyline curve.
   * \param q The point.
   * \return An index i such that q is in the x-range of cv[i].
   *         If q is not in the x-range of cv, returns(0xffffffff).
   */
  unsigned int locate_point(const X_monotone_curve_2 & cv,
                            const Point_2 & q) const
  {
    // First check whether the polyline curve really contains q in its x-range.
    unsigned int from = 0;
    Comparison_result res_from;
    unsigned int to = cv.size() - 1;
    Comparison_result res_to;

    Segment_traits_2 seg_traits;
    Compare_x_2 cmp_x = seg_traits.compare_x_2_object();
    Construct_min_vertex_2 min_ver = seg_traits.construct_min_vertex_2_object();
    Construct_max_vertex_2 max_ver = seg_traits.construct_max_vertex_2_object();

    res_from = cmp_x(min_ver(cv[from]), q);
    if (res_from == EQUAL) return from;
    
    res_to = cmp_x(max_ver(cv[to]), q);
    if (res_to == EQUAL) return to;
    
    if (res_from == res_to) return 0xffffffff;

    // Perform a binary search and locate the segment that contains q in its
    // x-range.
    unsigned int mid;
    Comparison_result res_mid_s, res_mid_t;

    while (to > from) {
      mid = (from + to) / 2;

      if (mid > from) {
	res_mid_s = cmp_x(min_ver(cv[mid]), q);

	if (res_mid_s == EQUAL) return mid;
	
	if (res_mid_s == res_from) from = mid;
	else to = mid - 1;
      } else {
	CGAL_assertion(mid < to);
	res_mid_t = cmp_x(max_ver(cv[mid]), q);

	if (res_mid_t == EQUAL) return mid;
	
	if (res_mid_t == res_to) to = mid;
	else from = mid + 1;
      }
    }

    // In case(from == to), and we know that the polyline contains the q:
    CGAL_assertion(from == to);
    return from;
  }

  /*! Find the index of the segment in the polyline that is defined to the
   * left(or to the right) of a the point q.
   * \param cv The polyline curve.
   * \param q The point.
   * \param to_right(true) if we wish to locate a segment to the right of q,
   *                (false) if we wish to locate a segment to its right.
   * \return An index i such that cv[i] is defined to the left(or to the
   *         right) of q, or(0xffffffff) if no such segment exists.
   */
  unsigned int locate_point_side(const X_monotone_curve_2 & cv,
                                 const Point_2 & q,
                                 const bool & to_right) const
  {
    // First locate a segment cv[i] that contains q in its x-range.
    Segment_traits_2 seg_traits;
    unsigned int i = locate_point(cv, q);

    if (i == 0xffffffff) return 0xffffffff;
   
    // If we seek an end-point to the right of q, q must be smaller than it.
    // If we seek an end-point to its left, q must be larger.
    const Comparison_result  cres = (to_right) ? SMALLER : LARGER;

    // Check whether x(q) is the x coordinate of an end-point of cv.
    if (seg_traits.compare_xy(q, seg_traits.curve_source(cv[i])) == EQUAL) {
      // If x(q) is at cv[i]'s source and its target is to the right(or left) 
      // of q, then cv[i] is defined to the right(or left) of q. 
      if (seg_traits.compare_xy(q, seg_traits.curve_target(cv[i])) == cres) {
	return i;
      }

      // Otherwise, since cv[i]'s source is cv[i-1]'s target, we should check
      // if the source of the previous curve is to the right(or left) of q.
      if (i > 0 &&
	  seg_traits.compare_x(q, seg_traits.curve_source(cv[i-1])) == cres)
      {
	return i-1;
      } else {
	return -1;
      }
    }
    else if (seg_traits.compare_xy(q, seg_traits.curve_target(cv[i])) == EQUAL)
    {
      // If x(q) is at cv[i]'s target and its source is to the right(or left) 
      // of q, then cv[i] is defined to the right(or left) of q. 
      if (seg_traits.compare_xy(q, seg_traits.curve_source(cv[i])) == cres)
      {
	return i;
      }

      // Otherwise, since cv[i]'s target is cv[i+1]'s source, we should check
      // if the target of the next curve is to the right(or left) of q.
      if (i <(cv.size() - 1) &&
	  seg_traits.compare_x(q, seg_traits.curve_target(cv[i+1])) == cres)
      {
	return i+1;
      } else {
	return 0xffffffff;
      }
    }

    // In case q is in cv[i]'s interior:
    return i;
  }

  /*! Check whether the cuve is defined from the left to the right.
   * \param cv The polyline curve.
   * \return (true) if the curve target is lexicographically larger than its
   *         source;(false) otherwise.
   */
  bool is_curve_to_right(const X_monotone_curve_2 & cv) const
  {
    const Point_2&    ps = seg_traits.curve_source(cv[0]);
    const Point_2&    pt = seg_traits.curve_target(cv[0]);
    Comparison_result res;
    
    res = seg_traits.compare_x(ps, pt);
    if (res == EQUAL)
      res = seg_traits.compare_xy(ps, pt);
				  
    return (res == SMALLER);
  }

  /*!
   * Find the nearest intersection point(or points) of two given curves to
   * the right(or to the left) lexicographically of a given point not 
   * including the point itself.
   * This function is used by both nearest_intersection_to_right() and
   * nearest_intersection_to_left() to avoid code duplication.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The refernece point.
   * \param to_right Should we search for an intersection to the right or to
   *                 the left of p.
   * \param p1 The first output point.
   * \param p2 The second output point.
   * \return 0 if cv1 and cv2 does not intersect to the right(or to the left,
   *           if to_right == false) of p;
   *         1 if they have a single intersection point, returned as p1;
   *         2 if they overlap, where p1, p2 are the endpoints of the overlap.
   */
  int _nearest_intersection_to_side(const X_monotone_curve_2 & cv1,
                                    const X_monotone_curve_2 & cv2,
                                    const Point_2 & p,
                                    const bool & to_right,
                                    Point_2 & p1, Point_2& p2) const
  {
    // Get the indices of the segments in cv1 and cv2 containing p.
    int i1 = locate_point(cv1, p);
    int i2 = locate_point(cv2, p);

    // Check if cv1 and cv2 are defined from left to right, and also
    // determine the desired comparison result(and its inverse).
    Comparison_result d_res;
    Comparison_result i_res;
    int n1 = cv1.size();
    int n2 = cv2.size();
    int inc1;
    int inc2;

    if (to_right) {
      d_res = SMALLER;
      i_res = LARGER;
      inc1 = is_curve_to_right(cv1) ? 1 : -1;
      inc2 = is_curve_to_right(cv2) ? 1 : -1;
    }
    else
    {
      d_res = LARGER;
      i_res = SMALLER;
      inc1 = is_curve_to_right(cv1) ? -1 : 1;
      inc2 = is_curve_to_right(cv2) ? -1 : 1;
    }

    // Check if the entire curve cv1 is defined to the right(left) of p. 
    if (i1 < 0) {
      if (inc1 == 1) {
	// The source is the leftmost(rightmost) point in cv1, so p must be 
	// to its left(right).
	// Otherwise, p is to the right(left) of the x-range of cv1. 
	if (seg_traits.compare_x(p, seg_traits.curve_source(cv1[0])) == d_res)
	  i1 = 0;
	else
	  return 0;
      } else {
        // if (inc1 == -1)
	// The target is the leftmost(rightmost) point in cv1, so p must be 
	// to its left(right).
	// Otherwise, p is to the right(left) of the x-range of cv1. 
	if (seg_traits.compare_x(p,
				 seg_traits.curve_target(cv1[n1-1])) == d_res)
	  i1 = n1-1;
	else
	  return 0;
      }
    }    

    // Check if the entire curve cv2 is defined to the right(left) of p. 
    if (i2 < 0) {
      if (inc2 == 1) {
	// The source is the leftmost(rightmost) point in cv2, so p must be 
	// to its left(right).
	// Otherwise, p is to the right(left) of the x-range of cv2. 
	if (seg_traits.compare_x(p, seg_traits.curve_source(cv2[0])) == d_res)
	  i2 = 0;
	else
	  return 0;
      } else {
        // if (inc2 == -1)
	// The target is the leftmost(rightmost) point in cv2, so p must be 
	// to its left(right).
	// Otherwise, p is to the right(left) of the x-range of cv2. 
	if (seg_traits.compare_x(p, seg_traits.curve_target(cv2[n2-1])) == d_res)
	  i2 = n2-1;
	else
	  return 0;
      }
    }

    // Try to locate the intersection point.
    Object obj;
    Segment_2 seg;
    Comparison_result res;
     
    while (i1 >= 0 && i1 < n1 && i2 >= 0 && i2 < n2) {      
      // Check if the two current segment intersect to the right(left) of p.
      if (to_right)
	obj = seg_traits.nearest_intersection_to_right(cv1[i1], cv2[i2], p);
      else
	obj = seg_traits.nearest_intersection_to_left(cv1[i1], cv2[i2], p);

      if (! obj.is_empty()) {
	// In case an overlap was detected, stop here:
	if (CGAL::assign(seg, obj)) {
	  p1 = seg_traits.curve_source(seg);
	  p2 = seg_traits.curve_target(seg);

	  return 2;
	}

        // The intersection is a single point:
        CGAL::assign(p1, obj);

	// In case we found a single intersection point, check whether it
	// is the next end-point of cv1[i1] or of cv2[i2].
	bool eq1, eq2;

	eq1 = seg_traits.point_equal(p1, 
				     (inc1 == 1) ? 
                                     seg_traits.curve_target(cv1[i1]) :
                                     seg_traits.curve_source(cv1[i1]));

	eq2 = seg_traits.point_equal(p1,
				     (inc2 == 1) ? 
                                     seg_traits.curve_target(cv2[i2]) :
                                     seg_traits.curve_source(cv2[i2]));
 	       
	if (!eq1 && !eq2) {
	  return 1;
	}
        
        // In case p1 equals one of the endpoints, simply assign this
        // endpoint to be p1.
        if (eq1) {
          p1 = (inc1 == 1) ? seg_traits.curve_target(cv1[i1]) :
            seg_traits.curve_source(cv1[i1]);
        }
        else { // if (eq2)
          p1 = (inc2 == 1) ? seg_traits.curve_target(cv2[i2]) :
            seg_traits.curve_source(cv2[i2]);
        }
        
	// Proceed to the next curves.
	if (eq1)
	  i1 += inc1;
	if (eq2)
	  i2 += inc2;

	if (i1 < 0 || i1 >= n1 || i2 < 0 || i2 >= n2)
	  return 1;

	// Check if the next curves overlap, and the nearest overlap endpoint
	// equals p1.
	Object   _obj;

	if (to_right)
	  _obj = seg_traits.nearest_intersection_to_right(cv1[i1], cv2[i2], p);
	else
	  _obj = seg_traits.nearest_intersection_to_left(cv1[i1], cv2[i2], p);

	if (CGAL::assign(seg, _obj))
	{	
	  Point_2 q1 = seg_traits.curve_source(seg);
	  Point_2 q2 = seg_traits.curve_target(seg);

	  if (seg_traits.point_equal(p1, q1)) {
	    // Now p1(== q1) --> q2 is an overlapping segment.
	    p2 = q2;
	    return 2;
	  } else if (seg_traits.point_equal(p1, q2)) {
	    // Now p1(== q2) --> q1 is an overlapping segment.
	    p2 = q1;
	    return 2;
	  }
	}

	// In this case we have a single intersection point.
	return 1;
      }

      // Find the segment whose end-point is the leftmost(rightmost) and move
      // forward(or backward) on its polyline.
      res = seg_traits.compare_x((inc1 == 1) ? 
                                 seg_traits.curve_target(cv1[i1]) :
                                 seg_traits.curve_source(cv1[i1]),
				 (inc2 == 1) ? 
                                 seg_traits.curve_target(cv2[i2]) :
                                 seg_traits.curve_source(cv2[i2]));

      if (res == d_res) {
	i1 += inc1;
      } else if (res == i_res) {
	i2 += inc2;
      } else {
	i1 += inc1;
	i2 += inc2;
      }
    }

    // No intersection found:
    return 0;
  }

};

/*!
 * A representation of a polyline, as used by the Arr_polyline_traits_2
 * traits class.
 */
template <class T_Segment_traits>
class Polyline_2 {
  friend class Arr_polyline_traits_2<T_Segment_traits>;

public:
  typedef T_Segment_traits                      Segment_traits_2;
  typedef typename Segment_traits_2::Point_2    Point_2;
  typedef typename Segment_traits_2::Curve_2    Segment_2;
  typedef ::std::vector<Segment_2>              Base;

private:
  // The segments that comprise the poyline:
  std::vector<Segment_2>                        segments;

public:
  /*! Default constructor. */
  Polyline_2() : segments() {}

  /*! Constructor from a range of points, defining the endpoints of the
   * polyline segments.
   * \param pts_begin An iterator pointing to the first point in the range.
   * \param pts_end An iterator pointing after the last point in the range.
   * \pre The are at least 2 points in the range.
   *      In other cases, an empty polyline will be created.
   */
  template <class InputIterator>
  Polyline_2(const InputIterator & pts_begin, const InputIterator & pts_end) :
    segments()
  {
    // Check if there are no points in the range:
    InputIterator  ps = pts_begin;

    if (ps == pts_end)
      return;

    // Construct a segment from each to adjacent points.
    InputIterator pt = ps;
    ++pt;

    while (pt != pts_end){
      segments.push_back(Segment_2(*ps, *pt));
      ++ps; ++pt;
    }
  }

  /* Append a point to the polyline */
  void push_back(const Point_2 & p)
  {
    Point_2 pt = p;
    Point_2 ps = (--segments.end()).source();
    segments.push_back(Segment_2(ps, pt));
  }

  /*! Create a bounding-box for the polyline.
   * \return The bounding-box.
   */
  Bbox_2 bbox() const
  {
    // Compute the union of the bounding boxes of all segments.
    unsigned int n = size();
    Bbox_2 bbox = (*this)[0].bbox();
    for (unsigned int i = 1; i < n; ++i)
      bbox = bbox +(*this)[i].bbox();
    return bbox;
  }

  class const_iterator;
  friend class const_iterator;

  /*!
   * An iterator for the polyline points.
   */
  class const_iterator
  {
  private:
    const Polyline_2<T_Segment_traits> *cvP;    // The polyline curve.
    int n_pts;                                  // Its number of points.
    T_Segment_traits seg_traits;                // Auxiliary variable.
    int i;                                      // The current point.
    bool is_forward;                            // Forward or reverse iterator.

    /*! Private constructor.
     * \param cv The scanned curve.
     * \param index The index of the segment.
     */
    const_iterator(const Polyline_2<T_Segment_traits>* _cvP,
                   const unsigned int & index,
                   const bool & _forward) :
      cvP(_cvP),
      i(index),
      is_forward(_forward)
    {
      if (cvP == NULL)
	n_pts = 0;
      else
	n_pts = (cvP->size() == 0) ? 0 :(cvP->size() + 1);
    }

    /*! ncrement the index.*/
    void increment()
    {
      if (cvP != NULL && i < n_pts) ++i;
    }

    /*! Decrement the index.*/
    void decrement()
    {
      if (cvP != NULL && i >= 0) --i;
    }

  public:
    /*! Default constructor. */
    const_iterator() :
      cvP(NULL),
      n_pts(0),
      i(0xffffffff),
      is_forward(true)
    {}

    /*! Dereference operator.
     * \return The current point.
     */
    Point_2 operator*() const
    {
      CGAL_assertion(cvP != NULL);
      CGAL_assertion(i >= 0 && i < n_pts);

      if (i == 0)
        // First point is the source of the first segment.
        return seg_traits.curve_source((*cvP)[0]);
      else
        // Return the target of the(i-1)'st segment
        return seg_traits.curve_target((*cvP)[i-1]);
    }

    /*! Increment operators. */
    void operator++() 
    {
      if (is_forward) increment();
      else decrement();
    }

    void operator++(int)
    {
      if (is_forward) increment();
      else decrement();
    }

    /*! Decrement operators. */
    void operator--()
    {
      if (is_forward) decrement();
      else increment();
    }

    void operator--(int)
    {
      if (is_forward) decrement();
      else increment();
    }

    /*! Equality operators. */
    bool operator==(const const_iterator& it) const
    {
      return (cvP == it.cvP && i == it.i);
    }

    bool operator!=(const const_iterator& it) const
    {
      return (cvP != it.cvP || i != it.i);
    }

    friend class Polyline_2<T_Segment_traits>;
  };

  /*! Get an iterator for the polyline points.
   * \return An iterator that points on the first point.
   */
  const_iterator begin() const
  {
    return (size() == 0) ? const_iterator(NULL, -1, true) :
      const_iterator(this, 0, true);
  }

  /*! Get a past-the-end iterator for the polyline points.
   * \return A past-the-end iterator.
   */
  const_iterator end() const
  {
    return (size() == 0) ? const_iterator(NULL, -1, true) :
      const_iterator(this, size() + 1, true);
  }

  /*! Get an reverse iterator for the polyline points.
   * \return An iterator that points on the last point.
   */
  const_iterator rbegin() const
  {
    return (size() == 0) ? const_iterator(NULL, -1, false) :
      const_iterator(this, size(), false);
  }

  /*! Get a reverse past-the-end iterator for the polyline points.
   * \return A reverse past-the-end iterator.
   */
  const_iterator rend() const
  {
    return (size() == 0) ? const_iterator(NULL, -1, false) :
      const_iterator(this, -1, false);
  }

  /*! Get the number of points contained in the polyline.
   * \return The number of points.
   */
  unsigned int points() const
  {
    return (size() == 0) ? 0 : size() + 1;
  }


private:

  /*! Append a segment to the polyline.
   * \param seg The new segment to be appended to the polyline.
   * \pre If the polyline is not empty, the segment source must be the
   * same as the target point of the last segment in the polyline.
   */
  inline void push_back(const Segment_2 & seg)
  {
    segments.push_back(seg);
  }

  /*! Get the number of segments that comprise the poyline.
   * \return The number of segments.
   */
  inline unsigned int size() const
  {
    return segments.size();
  }

  /*! Get the i'th segment of the polyline.
   * \param i The segment index(from 0 to size()-1).
   * \return A const reference to the segment.
   */
  inline const Segment_2 & operator[](const unsigned int & i) const
  {
    return segments[i];
  }

  /*! Clear the polyline. */
  inline void clear()
  {
    segments.clear();
  }
};

/*! Output operator for a polyline. */
template <class T_Segment_traits, class Stream_>
inline Stream_ & operator<<(Stream_ & os,
                            const Polyline_2<T_Segment_traits> & cv)
{
  typename Polyline_2<T_Segment_traits>::const_iterator ps = cv.begin();
  typename Polyline_2<T_Segment_traits>::const_iterator pt = ps; ++pt;

  while (pt != cv.end()) {
    typename T_Segment_traits::Curve_2 seg(*ps, *pt);
    os << seg;
    ++ps; ++pt;
  }
  return os;
}

/*! Specialized exporter for output stream.
 * In this case we export the number of points followed by the points
 */
template <class T_Segment_traits>
inline std::ostream & operator<<(std::ostream & os,
                                 const Polyline_2<T_Segment_traits> & pl)
{
  typedef Polyline_2<T_Segment_traits>  Curve_2;
  typename Curve_2::const_iterator it;

  // Print out the number of points in the polyline.
  os << pl.points();

  // Print out the polyline points.
  for (it = pl.begin(); it != pl.end(); ++it) 
    os << " " <<(*it);

  return os;
}

/*! Input operator for a polyline. */
template <class T_Segment_traits, class Stream_>
inline Stream_ & operator>>(Stream_ & is, Polyline_2<T_Segment_traits> & pl)
{
  typedef Polyline_2<T_Segment_traits>  Curve_2;
  typedef typename Curve_2::Point_2    Point_2;

  // Read the number of input points.
  int n_pts;

  is >> n_pts;

  // Read n_pts points to a list.
  Point_2 p;
  std::list<Point_2> pts;
  for (unsigned int i = 0; i < n_pts; ++i) {
    is >> p;
    pts.push_back(p);
  }

  // Create the polyline curve.
  pl = Curve_2(pts.begin(), pts.end());

  return is;
}

CGAL_END_NAMESPACE

#endif
