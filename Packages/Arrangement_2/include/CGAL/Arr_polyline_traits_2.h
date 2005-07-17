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
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>
//                 Ron Wein  <wein@post.tau.ac.il>

#ifndef CGAL_ARR_POLYLINE_TRAITS_2_H
#define CGAL_ARR_POLYLINE_TRAITS_2_H

/*! \file
 * The traits-class for the linear piece-wiese (polyline) type of curves of the
 * arrangement package.
 */

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_traits_2/Polyline_2.h>

CGAL_BEGIN_NAMESPACE

template <class SegmentTraits_>
class Arr_polyline_traits_2
{
public:

  typedef SegmentTraits_                                Segment_traits_2;

  // Tag defintion:
  typedef Tag_true                                      Has_left_category;
  typedef Tag_true                                      Has_merge_category;

protected:
  typedef Arr_polyline_traits_2<Segment_traits_2>       Self;

  // Data members:
  Segment_traits_2   seg_traits;           // The base segment-traits class.

public:

  /*! Default constructor */
  Arr_polyline_traits_2() :
    seg_traits()
  {}

  /// \name Types and functors inherited from the base segment traits.
  //@{

  // Traits types:
  typedef typename Segment_traits_2::Point_2            Point_2;
  typedef typename Segment_traits_2::Curve_2            Segment_2;

  typedef _Polyline_2<Segment_traits_2>                 Curve_2;
  typedef _X_monotone_polyline_2<Segment_traits_2>      X_monotone_curve_2;

  /*! Compare the x-coordinates of two points. */
  typedef typename Segment_traits_2::Compare_x_2        Compare_x_2;

  /*! Get a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object() const
  {
    return (seg_traits.compare_x_2_object());
  }

  /*! Compare two points lexigoraphically: by x, then by y. */
  typedef typename Segment_traits_2::Compare_xy_2       Compare_xy_2;
  
  /*! Get a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const
  {
    return (seg_traits.compare_xy_2_object());
  }

  ///@}

  /// \name Basic predicate functors (based on the segment traits).
  //@{

  class Construct_min_vertex_2
  {
  private:
    const Segment_traits_2  *seg_traits;

  public:

    /*! Constructor. */
    Construct_min_vertex_2 (const Segment_traits_2 *traits) :
      seg_traits (traits)
    {}

    /*!
     * Get the left endpoint of the x-monotone curve(segment).
     * \param cv The polyline curve.
     * \return The left endpoint.
     */
    const Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      CGAL_assertion (cv.size() > 0);

      return (seg_traits->construct_min_vertex_2_object()(cv[0]));
    }
  };
    
  /*! Get a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  {
    return Construct_min_vertex_2 (&seg_traits);
  }

  class Construct_max_vertex_2
  {
  private:
    const Segment_traits_2  *seg_traits;

  public:

    /*! Constructor. */
    Construct_max_vertex_2 (const Segment_traits_2 *traits) :
      seg_traits (traits)
    {}

    /*!
     * Get the right endpoint of the x-monotone curve(segment).
     * \param cv The polylinecurve.
     * \return The right endpoint.
     */
    const Point_2 operator() (const X_monotone_curve_2 & cv) const
    {
      CGAL_assertion(cv.size() > 0);

      return (seg_traits->construct_max_vertex_2_object() (cv[cv.size() - 1]));
    }
  };

  /*! Get a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  {
    return Construct_max_vertex_2 (&seg_traits);
  }

  class Is_vertical_2
  {
  private:
    const Segment_traits_2  *seg_traits;

  public:

    /*! Constructor. */
    Is_vertical_2 (const Segment_traits_2 *traits) :
      seg_traits (traits)
    {}

    /*!
     * Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2 & cv) const
    {
      // An x-monotone polyline can represent a vertical segment only if it
      // is comprised of a single vertical segment.
      return (cv.size() == 1 && seg_traits->is_vertical_2_object()(cv[0]));
    }
  };
  
  /*! Get an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const
  {
    return Is_vertical_2 (&seg_traits);
  }

  class Compare_y_at_x_2
  {
  private:
    const Segment_traits_2  *seg_traits;

  public:

    /*! Constructor. */
    Compare_y_at_x_2 (const Segment_traits_2 *traits) :
      seg_traits (traits)
    {}

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
      int     i = Self::_locate (seg_traits, cv, p);
      CGAL_precondition(i != -1);

      // Compare the segment cv[i] and p.
      return (seg_traits->compare_y_at_x_2_object() (p, cv[i]));
    }
  };
  friend class Compare_y_at_x_2;

  /*! Get a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  {
    return Compare_y_at_x_2 (&seg_traits);
  }

  class Compare_y_at_x_left_2
  {
  private:
    const Segment_traits_2  *seg_traits;

  public:

    /*! Constructor. */
    Compare_y_at_x_left_2 (const Segment_traits_2 *traits) :
      seg_traits (traits)
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
    Comparison_result operator() (const X_monotone_curve_2& cv1,
				  const X_monotone_curve_2& cv2,
				  const Point_2& p) const
    {
      // Get the indices of the segments in cv1 and cv2 containing p and
      // defined to its left.
      int     i1 = Self::_locate_side (seg_traits, cv1, p, false);
      int     i2 = Self::_locate_side (seg_traits, cv2, p, false);

      CGAL_precondition (i1 != -1);
      CGAL_precondition (i2 != -1);

      // Compare cv1[i1] and cv2[i2] at p's left.
      return (seg_traits->compare_y_at_x_left_2_object()(cv1[i1], cv2[i2], p));
    }
  };
  friend class Compare_y_at_x_left_2;

  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  {
    return Compare_y_at_x_left_2 (&seg_traits);
  }

  class Compare_y_at_x_right_2
  {
  private:
    const Segment_traits_2  *seg_traits;

  public:

    /*! Constructor. */
    Compare_y_at_x_right_2 (const Segment_traits_2 *traits) :
      seg_traits (traits)
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
    Comparison_result operator() (const X_monotone_curve_2 & cv1,
				  const X_monotone_curve_2 & cv2,
				  const Point_2 & p) const
    {
      // Get the indices of the segments in cv1 and cv2 containing p and
      // defined to its right.
      int     i1 = Self::_locate_side (seg_traits, cv1, p, true);
      int     i2 = Self::_locate_side (seg_traits, cv2, p, true);

      CGAL_precondition (i1 != -1);
      CGAL_precondition (i2 != -1);

      // Compare cv1[i1] and cv2[i2] at p's right.
      return (seg_traits->compare_y_at_x_right_2_object()(cv1[i1],cv2[i2], p));
    }
  };
  friend class Compare_y_at_x_right_2;

  /*! Get a Compare_y_at_x_right_2 functor object. */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  {
    return Compare_y_at_x_right_2 (&seg_traits);
  }

  class Equal_2
  {
  private:
    const Segment_traits_2  *seg_traits;

  public:

    /*! Constructor. */
    Equal_2 (const Segment_traits_2 *traits) :
      seg_traits (traits)
    {}

    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same;(false) otherwise.
     */
    bool operator() (const Point_2& p1, const Point_2& p2) const
    {
      return (seg_traits->equal_2_object() (p1, p2));
    }

    /*!
     * Check if the two x-monotone curves are the same(have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same;(false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
		     const X_monotone_curve_2& cv2) const
    {
      // The two curves must contain the same number of segments.
      unsigned int n1 = cv1.size();
      unsigned int n2 = cv2.size();
      if (n1 != n2)
	return (false);

      // Check the pairwise equality of the contained segments.
      typename Segment_traits_2::Equal_2 equal = seg_traits->equal_2_object();
      unsigned int                       i;

      for (i = 0; i < n1; i++)
      {
        if (!equal(cv1[i], cv2[i]))
	  return (false);
      }
      return (true);
    }
  };

  /*! Get an Equal_2 functor object. */
  Equal_2 equal_2_object() const
  {
    return Equal_2 (&seg_traits);
  }
  ///@}

  /// \name Construction functors (based on the segment traits).
  //@{

  class Make_x_monotone_2
  {
  private:
    const Segment_traits_2  *seg_traits;

  public:

    /*! Constructor. */
    Make_x_monotone_2 (const Segment_traits_2 *traits) :
      seg_traits (traits)
    {}

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
    OutputIterator operator() (const Curve_2& cv, OutputIterator oi) const
    { 
      // Go over all points in the input curve.
      typename Curve_2::const_iterator       ps = cv.begin();
      typename Curve_2::const_iterator       end = cv.end();

      if (ps == end)
	// Empty polyline:
	return (oi);

      typename Curve_2::const_iterator       pt = ps;
      ++pt;

      if (pt == end)
      {
	// The polyline contains a single isolated point:
	*oi = make_object (*ps);
	++oi;
	return (oi);
      }

      // Locate points where the x-order changes:
      typename Segment_traits_2::Compare_x_2     compare_x =
        seg_traits->compare_x_2_object();
      typename Segment_traits_2::Compare_xy_2    compare_xy =
        seg_traits->compare_xy_2_object();
      Comparison_result                          x_res;
      Comparison_result                          xy_res;
      typename Curve_2::const_iterator           x_begin = ps;
      Comparison_result                          curr_x_res;
      Comparison_result                          curr_xy_res;
     
      x_res = compare_x (*ps, *pt);
      if (x_res != EQUAL)
	xy_res = x_res;
      else
	xy_res = compare_xy (*ps, *pt);

      ++ps; ++pt;
      while (pt != end)
      {
	curr_x_res = compare_x (*ps, *pt);
	if (curr_x_res != EQUAL)
	  curr_xy_res = curr_x_res;
	else
	  curr_xy_res = compare_xy (*ps, *pt);

	if (curr_x_res != x_res || curr_xy_res != xy_res)
	{
	  // Create a new x-monotone polyline from the range of points
	  // [x_begin, pt):
	  *oi = make_object (X_monotone_curve_2 (x_begin, pt));
	  ++oi;

	  x_begin = ps;
	  x_res = curr_x_res;
	  xy_res = curr_xy_res;
	}

	++ps; ++pt;
      }

      // Create an x-monotone polyline from the remaining points.
      CGAL_assertion (x_begin != end);
      *oi = make_object (X_monotone_curve_2 (x_begin, end));
      ++oi;
      
      return (oi);
    }
  };
  
  /*! Get a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  {
    return Make_x_monotone_2 (&seg_traits);
  }

  class Split_2
  {
  private:
    Segment_traits_2        *seg_traits;

  public:

    /*! Constructor. */
    Split_2 (Segment_traits_2 *traits) :
      seg_traits (traits)
    {}

  public:

    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve(p is its right endpoint).
     * \param c2 Output: The right resulting subcurve(p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator() (const X_monotone_curve_2& cv,
		     const Point_2& p,
		     X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
        seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Equal_2                equal =
	seg_traits->equal_2_object();

      // Make sure the split point is not one of the curve endpoints.
      CGAL_precondition (!equal(min_vertex(cv[0]), p));
      CGAL_precondition (!equal(max_vertex(cv[cv.size()-1]), p));

      // Locate the segment on the polyline cv that contains p.
      int     i = Self::_locate (seg_traits, cv, p);
      CGAL_precondition (i != -1);

      // Clear the output curves.
      c1.clear(); 
      c2.clear();

      // Push all segments labeled(0, 1, ... , i-1) into c1.
      int     j;
      for (j = 0; j < i; j++)
	c1.push_back(cv[j]);

      // Check whether the split point is cv[i]'s source of target.
      if (equal(max_vertex(cv[i]), p))
      {
	// The entire i'th segment belongs to c1:
	c1.push_back (cv[i]);
      }
      else if (equal(min_vertex(cv[i]), p))
      {
	// The entire i'th segments belongs to c2:
	c2.push_back (cv[i]);
      }
      else
      {
        // The i'th segment should be split: The left part (seg1) goes to cv1,
	// and the right part (seg2) goes to cv2.
        Segment_2   seg1, seg2;
        seg_traits->split_2_object() (cv[i], p, seg1, seg2);

        c1.push_back(seg1);
        c2.push_back(seg2);
      }

      // Push all segments labeled(i+1, i+2, ... , n-1) into cv1.
      int    n = static_cast<int> (cv.size());

      for (j = i+1; j < n; j++)
	c2.push_back(cv[j]);

      return;
    }
  };
  friend class Split_2;

  /*! Get a Split_2 functor object. */
  Split_2 split_2_object()
  {
    return Split_2 (&seg_traits);
  }

  class Intersect_2
  {
  private:
    Segment_traits_2        *seg_traits;

  public:

    /*! Constructor. */
    Intersect_2 (Segment_traits_2 *traits) :
      seg_traits (traits)
    {}

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
    OutputIterator operator() (const X_monotone_curve_2& cv1,
			       const X_monotone_curve_2& cv2,
			       OutputIterator oi)
    {
      typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
        seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Equal_2                equal =
	seg_traits->equal_2_object();
      typename Segment_traits_2::Compare_xy_2           compare_xy =
        seg_traits->compare_xy_2_object();
      typename Segment_traits_2::Intersect_2            intersect =
        seg_traits->intersect_2_object();
      typename Segment_traits_2::Compare_y_at_x_2       compare_y_at_x =
        seg_traits->compare_y_at_x_2_object();
      
      const int           n1 = static_cast<int> (cv1.size());
      const int           n2 = static_cast<int> (cv2.size());
      int                 i1 = 0;
      int                 i2 = 0;
      X_monotone_curve_2  ocv;           // Used to represent overlaps.

      Comparison_result   left_res = compare_xy (min_vertex(cv1[i1]),
						 min_vertex(cv2[i2]));
      
      if (left_res == SMALLER)
      {
	// cv1's left endpoint is to the left of cv2's left endpoint:
	// Locate the index i1 of the segment in cv1 which contains cv2's
	// left endpoint.
        i1 = Self::_locate (seg_traits, cv1, min_vertex(cv2[i2]));

	if (i1 == -1)
	  return (oi);

        if (equal (max_vertex(cv1[i1]), min_vertex(cv2[i2])))
	{
          i1++;
          left_res = EQUAL;
        }
      }
      else if (left_res == LARGER)
      {
	// cv1's left endpoint is to the right of cv2's left endpoint:
	// Locate the index i2 of the segment in cv2 which contains cv1's
	// left endpoint.
	i2 = Self::_locate (seg_traits, cv2, min_vertex(cv1[i1]));

	if (i2 == -1)
	  return (oi);

        if (equal (max_vertex(cv2[i2]), min_vertex(cv1[i1])))
	{
          i2++;
          left_res = EQUAL;
        }
      }

      // Check if the the left endpoint lies on the other polyline.
      bool               left_coincides = (left_res == EQUAL);
      bool               left_overlap = false;

      if (left_res == SMALLER)
      {
        left_coincides = (compare_y_at_x (min_vertex(cv2[i2]),
					  cv1[i1]) == EQUAL);
      }
      else if (left_res == LARGER)
      {
        left_coincides = (compare_y_at_x (min_vertex(cv1[i1]),
					  cv2[i2]) == EQUAL);
      }

      // The main loop: Go simultaneously over both polylines.
      Comparison_result  right_res;
      bool               right_coincides;
      bool               right_overlap;

      while (i1 < n1 && i2 < n2)
      {
        right_res = compare_xy (max_vertex(cv1[i1]), max_vertex(cv2[i2]));

        right_coincides = (right_res == EQUAL);
        if (right_res == SMALLER)
	{
          right_coincides = (compare_y_at_x (max_vertex(cv1[i1]),
					     cv2[i2]) == EQUAL);
	}
        else if (right_res == LARGER)
	{
          right_coincides = (compare_y_at_x (max_vertex(cv2[i2]),
					     cv1[i1]) == EQUAL);
	}

        right_overlap = false;
        
	if (!right_coincides && !left_coincides)
	{
          // Non of the endpoints of the current segment of one polyline
          // coincides with the curent segment of the other polyline:
          // Output the intersection if exists.
          oi = intersect(cv1[i1], cv2[i2], oi);
        }
	else if (right_coincides && left_coincides)
	{
          // An overlap exists between the current segments of the polylines:
          // Output the overlapping segment.
          right_overlap = true;
          if (left_res == SMALLER)
	  {
            if (right_res == SMALLER)
	    {
              Segment_2 seg (min_vertex(cv2[i2]), max_vertex(cv1[i1]));
              ocv.push_back (seg);
            }
	    else
	    {
              Segment_2 seg (min_vertex(cv2[i2]), max_vertex(cv2[i2]));
              ocv.push_back (seg);
            }
          }
	  else
	  {
           if (right_res == SMALLER)
	   {
              Segment_2 seg (min_vertex(cv1[i1]), max_vertex(cv1[i1]));
              ocv.push_back (seg);
	   }
	   else
	   {
              Segment_2 seg (min_vertex(cv1[i1]), max_vertex(cv2[i2]));
              ocv.push_back (seg);
            }
          }
        }
	else if (left_coincides && !right_coincides)
	{
          // The left point of the current segment of one polyline
          // coincides with the current segment of the other polyline.
          if (left_overlap)
	  {
            // An overlap occured at the previous iteration:
            // Output the overlapping polyline.
            CGAL_assertion(ocv.size() > 0);
            *oi = make_object(ocv);
            ++oi;
	    ocv.clear();
          }
	  else
	  {
            // The left point of the current segment of one polyline
            // coincides with the current segment of the other polyline, and
            // no overlap occured at the previous iteration:
            // Output the intersection point. The derivative of at least one of
	    // the polylines is not defined at this point, so we give it
	    // multiplicity 0.
            if (left_res == SMALLER)
	    {
              std::pair<Point_2, unsigned int>  p (min_vertex(cv2[i2]), 0);
              *oi = make_object (p);
	      ++oi;
            }
	    else
	    {
              std::pair<Point_2, unsigned int>  p (min_vertex(cv1[i1]), 0);
              *oi = make_object (p);
	      ++oi;
            }
          }
        }

	// Proceed forward.
        if (right_res != SMALLER)
	  ++i2;
        if (right_res != LARGER)
	  ++i1;

        if (right_res == SMALLER)
	  left_res = LARGER;
        else if (right_res == LARGER)
	  left_res = SMALLER;

        left_coincides = right_coincides;
        left_overlap = right_overlap;
      }

      // Output the remaining overlapping polyline, if necessary.
      if (ocv.size() > 0)
      {
        *oi = make_object(ocv);
        ++oi;
      }

      // RWRW: This is a patch, for Efi to fix.
      // Somehow, when the two right endpoints are equal, we do not get
      // an intersection point back.
      if (compare_xy (max_vertex(cv1[n1 - 1]),
		      max_vertex(cv2[n2 - 1])) == EQUAL)
      {
	std::pair<Point_2, unsigned int> ip (max_vertex(cv1[n1 - 1]), 0);
	*oi = make_object (ip);
	++oi;
      }

      return (oi);
    }
  };
  friend class Intersect_2;
  
  /*! Get an Intersect_2 functor object. */
  Intersect_2 intersect_2_object()
  {
    return Intersect_2 (&seg_traits);
  }

  class Are_mergeable_2
  {
  private:
    const Segment_traits_2  *seg_traits;

  public:

    /*! Constructor. */
    Are_mergeable_2 (const Segment_traits_2 *traits) :
      seg_traits (traits)
    {}

    /*!
     * Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable, that is, they share a
     * common endpoint; (false) otherwise.
     */
    bool operator()(const X_monotone_curve_2 & cv1,
                    const X_monotone_curve_2 & cv2) const
    {
      typename Segment_traits_2::Construct_min_vertex_2    min_vertex =
        seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Equal_2                   equal =
	seg_traits->equal_2_object();
      
      const unsigned int  n1 = cv1.size();
      const unsigned int  n2 = cv2.size();
      
      return (equal (max_vertex(cv1[n1 - 1]), min_vertex(cv2[0])) ||
	      equal (max_vertex(cv2[n2 - 1]), min_vertex(cv1[0])));
    }
  };
  
  /*! Get an Are_mergeable_2 functor object. */
  Are_mergeable_2 are_mergeable_2_object() const
  {
    return Are_mergeable_2 (&seg_traits);
  }

  class Merge_2
  {
  private:
    Segment_traits_2        *seg_traits;

  public:

    /*! Constructor. */
    Merge_2 (Segment_traits_2 *traits) :
      seg_traits (traits)
    {}

  public:

    /*!
     * Merge two given x-monotone curves into a single curve(segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable. That is, they share a common
     *      endpoint.
     */
    void operator() (const X_monotone_curve_2& cv1,
		     const X_monotone_curve_2& cv2,
		     X_monotone_curve_2& c) const
    {
      typename Segment_traits_2::Construct_min_vertex_2    min_vertex =
        seg_traits->construct_min_vertex_2_object();
      typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
        seg_traits->construct_max_vertex_2_object();
      typename Segment_traits_2::Equal_2                   equal =
	seg_traits->equal_2_object();
      
      const unsigned int  n1 = cv1.size();
      const unsigned int  n2 = cv2.size();
      unsigned int        i;

      c.clear();
      if (equal (max_vertex(cv1[n1 - 1]), min_vertex(cv2[0])))
      {
	// cv2 extends cv1 to the right:
	for (i = 0; i < n1 - 1; ++i)
	  c.push_back(cv1[i]);

	// Try to merge tthe to contiguous line segments:
	if (seg_traits->are_mergeable_2_object() (cv1[n1 - 1], cv2[0]))
	{
	  Segment_2       seg;
	  seg_traits->merge_2_object() (cv1[n1 - 1], cv2[0], seg);
	  c.push_back (seg);
	}
	else
	{
	  c.push_back (cv1[n1 - 1]);
	  c.push_back (cv2[0]);
	}

	for (i = 1; i < n2; ++i)
	  c.push_back(cv2[i]);
      }
      else if (max_vertex(cv2[n2 - 1]), min_vertex(cv1[0]))
      {
	// cv1 extends cv2 to the right:
	for (i = 0; i < n2 - 1; ++i)
	  c.push_back(cv2[i]);

	// Try to merge tthe to contiguous line segments:
	if (seg_traits->are_mergeable_2_object() (cv2[n2 - 1], cv1[0]))
	{
	  Segment_2       seg;
	  seg_traits->merge_2_object() (cv2[n2 - 1], cv1[0], seg);
	  c.push_back (seg);
	}
	else
	{
	  c.push_back (cv2[n2 - 1]);
	  c.push_back (cv1[0]);
	}

	for (i = 1; i < n1; ++i)
	  c.push_back(cv1[i]);
      }
      else
      {
	CGAL_precondition_msg (false, "The curves are not mergeable.");
      }

      return;
    }
  };
  
  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object() const
  {
    return Merge_2 (&seg_traits);
  }
  ///@}
  
  /// \name functors required for the landmark point location

  //typedef typename Segment_traits_2::Approximate_number_type
  //  Approximate_number_type;
  //typedef typename Segment_traits_2::Approximate_2            Approximate_2;
  
  /*! Get an Approximate_2 functor object. */
  //Approximate_2 approximate_2_object() const
  //{
  //  return Approximate_2();
  //}
  ///@}

private:
  
  /*!
   * Return the index of the segment in the polyline that contains the
   * point q in its x-range. The function performs a binary search, so if the
   * point q is in the x-range of the polyline with n segments, the segment
   * containing it can be located in O(log n) operations.
   * \param cv The polyline curve.
   * \param q The point.
   * \return An index i such that q is in the x-range of cv[i].
   *         If q is not in the x-range of cv, returns(-1).
   */
  static int _locate (const Segment_traits_2 *p_seg_traits,
		      const X_monotone_curve_2& cv,
		      const Point_2& q)
  {
    // First check whether the polyline curve really contains q in its x-range.
    int               from = 0;
    Comparison_result res_from;
    int               to = static_cast<int>(cv.size()) - 1;
    Comparison_result res_to;

    typename Segment_traits_2::Compare_x_2            compare_x =
      p_seg_traits->compare_x_2_object();
    typename Segment_traits_2::Construct_min_vertex_2 min_vertex =
      p_seg_traits->construct_min_vertex_2_object();
    typename Segment_traits_2::Construct_max_vertex_2 max_vertex =
      p_seg_traits->construct_max_vertex_2_object();

    res_from = compare_x (min_vertex (cv[from]), q);
    if (res_from == EQUAL)
      return (from);
    
    res_to = compare_x (max_vertex (cv[to]), q);
    if (res_to == EQUAL)
      return (to);
    
    if (res_from == res_to)
      // q is not in the x-range of cv:
      return (-1);

    // Perform a binary search and locate the segment that contains q in its
    // x-range.
    int                mid;
    Comparison_result  res_mid_s, res_mid_t;

    while (to > from)
    {
      mid = (from + to) / 2;

      if (mid > from)
      {
	res_mid_s = compare_x (min_vertex (cv[mid]), q);

	if (res_mid_s == EQUAL)
	  return (mid);
	
	if (res_mid_s == res_from)
	  from = mid;
	else
	  to = mid - 1;
      } 
      else
      {
	CGAL_assertion(mid < to);
	res_mid_t = compare_x (max_vertex (cv[mid]), q);

	if (res_mid_t == EQUAL)
	  return (mid);
	
	if (res_mid_t == res_to)
	  to = mid;
	else
	  from = mid + 1;
      }
    }

    // In case (from == to), and we know that the polyline contains the q:
    CGAL_assertion(from == to);
    return (from);
  }

  /*!
   * Find the index of the segment in the polyline that is defined to the
   * left (or to the right) of the point q.
   * \param cv The polyline curve.
   * \param q The point.
   * \param to_right (true) if we wish to locate a segment to the right of q,
   *                (false) if we wish to locate a segment to its right.
   * \return An index i such that segments[i] is defined to the left (or to the
   *         right) of q, or(-1) if no such segment exists.
   */
  static int _locate_side (const Segment_traits_2 *p_seg_traits,
			   const X_monotone_curve_2& cv,
                           const Point_2& q, const bool& to_right)
  {
    // First locate a segment segments[i] that contains q in its x-range.
    int     i = _locate (p_seg_traits, cv, q);
    if (i == -1)
      return (-1);
   
    typename Segment_traits_2::Equal_2 equal = p_seg_traits->equal_2_object();

    if (equal (p_seg_traits->construct_min_vertex_2_object() (cv[i]), q))
    {
      // q is the left endpoint of the i'th segment:
      if (to_right)
	return (i);
      else if (i == 0)
	return (-1);
      else
	return (i - 1);
    }

    if (equal (p_seg_traits->construct_max_vertex_2_object() (cv[i]), q))
    {
      // q is the right endpoint of the i'th segment:
      if (!to_right)
	return (i);
      else if (i == static_cast<int>(cv.size()) - 1)
	return (-1);
      else
	return (i + 1);
    }    
    
    // In case q is in cv[i]'s interior:
    return (i);
  }
};


CGAL_END_NAMESPACE

#endif
