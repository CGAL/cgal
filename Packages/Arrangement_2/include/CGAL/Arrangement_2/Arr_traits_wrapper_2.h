// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein             <wein@post.tau.ac.il>
//                 (based on old version by Iddo Hanniel
//                                          Eyal Flato
//                                          Oren Nechushtan
//                                          Efi Fogel
//                                          Ron Wein
//                                          Idit Haran)
#ifndef CGAL_ARR_TRAITS_WRAPPER_2_H
#define CGAL_ARR_TRAITS_WRAPPER_2_H

/*! \file
 * Definitions of the wrapper classes for the arrangement traits class.
 */

#include <CGAL/config.h>
#include <CGAL/tags.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A traits-class wrapper that extends the basic traits-class interface.
 */
template <class ArrangementBasicTraits_>
class Arr_traits_basic_wrapper_2 : public ArrangementBasicTraits_
{
public:

  // Traits-class geometric types.
  typedef ArrangementBasicTraits_               Base;
  typedef typename Base::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Base::Point_2                Point_2;

  // Tags.
  typedef typename Base::Has_left_category      Has_left_category;

  /// \name Construction.
  //@{
  /*! Default constructor. */
  Arr_traits_basic_wrapper_2 () :
    Base()
  {}

  /*! Constructor from a base-traits class. */
  Arr_traits_basic_wrapper_2 (const Base& traits) :
    Base (traits)
  {}
  //@}

  // Inherited functors:
  typedef typename Base::Compare_x_2            Compare_x_2;
  typedef typename Base::Compare_xy_2           Compare_xy_2;
  typedef typename Base::Construct_min_vertex_2 Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2 Construct_max_vertex_2;
  typedef typename Base::Is_vertical_2          Is_vertical_2;
  typedef typename Base::Compare_y_at_x_2       Compare_y_at_x_2;
  typedef typename Base::Compare_y_at_x_right_2 Compare_y_at_x_right_2;
  typedef typename Base::Equal_2                Equal_2;

  /// \name Overriden functors.
  //@{
  class Compare_y_at_x_left_2
  {
  public:
    /*!
     * Compare two curves immediately to the left of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The query point.
     * \pre The two curves intersect at p, and they are defined to its left.
     * \return SMALLER if cv1 lies below cv2 to the left of q;
     *         LARGER if cv1 lies above cv2 to the left of q;
     *         EQUAL in case of an overlap to the left of q.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p) const 
    {
      // The function is implemented based on the Has_left category. If the 
      // category indicates that the "left" version is available, it calls the
      // function with same name defined in the base class. Otherwise, it 
      // uses other predicates to provide this comparison.
      return _compare_y_at_x_left_imp (cv1, cv2, p, Has_left_category());
    }

  private:

    /*!
     * Implementation of the operator() in case the HasLeft tag is true.
     */
    Comparison_result _compare_y_at_x_left_imp (const X_monotone_curve_2& cv1,
                                                const X_monotone_curve_2& cv2,
                                                const Point_2& p,
                                                Tag_true) const
    {
      Base                    tr;
      return (tr.compare_y_at_x_left_2_object() (cv1, cv2, p));
    }

    /*!
     * Implementation of the operator() in case the HasLeft tag is false.
     */
    Comparison_result _compare_y_at_x_left_imp (const X_monotone_curve_2& cv1,
                                                const X_monotone_curve_2& cv2,
                                                const Point_2& p,
                                                Tag_false) const
    {
      Base                    tr;
      Construct_min_vertex_2  min_vertex = tr.construct_min_vertex_2_object();
      Equal_2                 equal = tr.equal_2_object();

      // Make sure that p lies on both curves, and that both are defined to its
      // right (so their right endpoint is lexicographically larger than p).
      CGAL_precondition_code (
        Compare_xy_2          compare_xy = tr.compare_xy_2_object();
        Compare_y_at_x_2      compare_y_at_x = tr.compare_y_at_x_2_object();
      );

      CGAL_precondition (compare_y_at_x (p, cv1) == EQUAL &&
                         compare_y_at_x (p, cv2) == EQUAL);

      CGAL_precondition (compare_xy(min_vertex (cv1), p) == SMALLER &&
                         compare_xy(min_vertex (cv2), p) == SMALLER);

      // If one of the curves is vertical, it is below the other one.
      Is_vertical_2           is_vertical = tr.is_vertical_2_object();

      if (is_vertical(cv1))
      {
        if (is_vertical (cv2))
          // Both are vertical:
          return (EQUAL);
        else
          return (SMALLER);
      }
      else if (is_vertical (cv2))
      {
        return (LARGER);
      }

      // Get the left endpoints of cv1 and cv2.
      Point_2        left1 = min_vertex(cv1);
      Point_2        left2 = min_vertex(cv2);

      if (equal (left1, left2))
      {
        // The two curves have a common left endpoint:
        // Compare them to the right of this point.
        return (tr.compare_y_at_x_right_2_object() (cv1, cv2, left1));
      }    

      // Compare the relative position of the curves at the righmost of left1
      // and left2:
      Compare_y_position_2   compare_y_position;
      return (compare_y_position (cv1, cv2));
    }
  };

  /*! Get a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object () const
  {
    return Compare_y_at_x_left_2();
  }
  //@}

  /// \name Additional auxiliary functors.
  //@{
  class Is_in_x_range_2
  {
  public:
    /*!
     * Check whether the given point is in the x-range of the given x-monotone
     * curves.
     * \param cv The x-monotone curve.
     * \param p The point.
     * \return (true) if x(cv_left) <= x(p) <= x(cv_right); (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv, const Point_2& p) const
    {
      // Compare p to the left and right endpoints of the curve.
      Base                    tr;
      Construct_min_vertex_2  min_vertex = tr.construct_min_vertex_2_object();
      Construct_max_vertex_2  max_vertex = tr.construct_max_vertex_2_object();
      Compare_x_2             compare_x = tr.compare_x_2_object();

      if (compare_x (p, min_vertex(cv)) == SMALLER)
        // p is to the left of the x-range.
        return (false);

      if (compare_x (p, max_vertex(cv)) == LARGER)
        // p is to the right of the x-range.
        return (false);

      // p is in the x-range of cv:
      return (true);
    }

    /*!
     * Check whether the x-ranges of the given x-monotone curves overlap.
     * \param cv1 The first x-monotone curve.
     * \param cv2 The second x-monotone curve.
     * \return (true) if there is an overlap in the x-ranges of the given
     *         curves; (false) otherwise.
     */
    bool operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2) const
    {
      // Find p_l, the rightmost of the two left endpoints of the curves
      // and p_r, the leftmost of the two right endpoints of the curves.
      Base                    tr;
      Construct_min_vertex_2  min_vertex = tr.construct_min_vertex_2_object();
      Construct_max_vertex_2  max_vertex = tr.construct_max_vertex_2_object();
      Compare_x_2             compare_x = tr.compare_x_2_object();

      const Point_2&  l1 = min_vertex(cv1);
      const Point_2&  r1 = max_vertex(cv1);
      const Point_2&  l2 = min_vertex(cv2);
      const Point_2&  r2 = max_vertex(cv2);
      const Point_2&  p_l = ((compare_x (l1, l2) == LARGER) ? l1 : l2);
      const Point_2&  p_r = ((compare_x (r1, r2) == SMALLER) ? r1 : r2);

      // The two curves overlap in their x-range if and only if p_l is not
      // to the right of p_r.
      return (compare_x (p_l, p_r) != LARGER);
    }
  };

  /*! Get an Is_in_x_range_2 functor object. */
  Is_in_x_range_2 is_in_x_range_2_object () const
  {
    return Is_in_x_range_2();
  }

  class Compare_y_position_2
  {
  public:
    /*!
     * Get the relative of two x-monotone curves with overlapping x-ranges
     * that are disjoint in their interiors.
     * \param cv1 The first x-monotone curve.
     * \param cv2 The second x-monotone curve.
     * \pre The x-ranges of the two curves overlap.
     * \return SMALLER if cv1 lies below cv2;
     *         LARGER if cv1 lies above cv2;
     *         EQUAL in case of an overlap (illegal input).
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2) const
    {
      CGAL_precondition_code (
        Is_in_x_range_2  is_in_x_range;
      );
      CGAL_precondition (is_in_x_range (cv1, cv2));

      // Get the left endpoints of cv1 and cv2.
      Base                    tr;
      Construct_min_vertex_2  min_vertex = tr.construct_min_vertex_2_object();
      Construct_max_vertex_2  max_vertex = tr.construct_max_vertex_2_object();
      Compare_xy_2            compare_xy = tr.compare_xy_2_object();
      Compare_y_at_x_2        compare_y_at_x = tr.compare_y_at_x_2_object();
      Compare_y_at_x_right_2  compare_y_at_x_right =
                                           tr.compare_y_at_x_right_2_object();

      const Point_2&  left1 = min_vertex(cv1);
      const Point_2&  left2 = min_vertex(cv2);

      // Locate the rightmost point of left1 and left2 and compare its position
      // to the other curve.
      Comparison_result  l_res = compare_xy (left1, left2);
      Comparison_result  res;

      if (l_res != SMALLER)
      {
        // left1 is in the x-range of cv2:
        res = compare_y_at_x (left1, cv2);

        if (res == EQUAL)
        {
          // The two curves intersect at left1. If both curves are defined to
          // the right of the reference point, we can compare them to its 
          // right. Otherwise, their share a common endpoint (which is the only
          // overlap in their x-ranges) and are really equal.
          if (l_res == EQUAL)
            res = compare_y_at_x_right (cv1, cv2, left1);
        }

        return (res);
      }
      else
      {
        // left2 is in the x-range of cv1:
        res = compare_y_at_x (left2, cv1);

        if (res == EQUAL)
        {
          // The two curves share a common endpoint (which is the only overlap
          // in their x-ranges) and are really equal.
          return (EQUAL);
        }

        // Swap the result:
        return ((res == SMALLER) ? LARGER : SMALLER);
      }
    }
  };

  /*! Get a Compare_y_position_2 functor object. */
  Compare_y_position_2 compare_y_position_2_object () const
  {
    return Compare_y_position_2();
  }

  class Is_between_cw_2
  {
  public:
    /*!
     * Check if the given query curve is encountered when rotating the first
     * curve in a clockwise direction around a given point until reaching the
     * second curve.
     * \param cv The query curve.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The point around which we rotate cv1.
     * \param cv_equal_cv1 Output: does cv equal cv1.
     * \param cv_equal_cv2 Output: does cv equal cv2.
     * \pre p is an end-point of all three curves.
     * \return (true) if cv is between cv1 and cv2; (false) otherwise.
     *         If cv overlaps cv1 or cv2 the result is always (false).
     *         If cv1 and cv2 overlap, the result is (true), unless cv 
     *         also overlaps them.
     */
    bool operator() (const X_monotone_curve_2& cv, 
                     const X_monotone_curve_2& cv1, 
                     const X_monotone_curve_2& cv2, 
                     const Point_2& p,
                     bool& cv_equal_cv1, 
                     bool& cv_equal_cv2) const
    {
      Base                    tr;
      Construct_min_vertex_2  min_vertex = tr.construct_min_vertex_2_object();
      Construct_max_vertex_2  max_vertex = tr.construct_max_vertex_2_object();
      Equal_2                 equal = tr.equal_2_object();
      Compare_y_at_x_right_2  compare_y_at_x_right =
                                           tr.compare_y_at_x_right_2_object();
      Compare_y_at_x_left_2   compare_y_at_x_left;

      // Find the direction of each $x$-monotone curve with respect to p.
      bool               cv_to_right;
      bool               cv1_to_right;
      bool               cv2_to_right;

      if (equal (min_vertex (cv), p))
      {
        cv_to_right = true;
      }
      else
      {
        CGAL_precondition (equal (max_vertex (cv), p));
        cv_to_right = false;        
      }

      if (equal (min_vertex (cv1), p))
      {
        cv1_to_right = true;
      }
      else
      {
        CGAL_precondition (equal (max_vertex (cv1), p));
        cv1_to_right = false;        
      }

      if (equal (min_vertex (cv2), p))
      {
        cv2_to_right = true;
      }
      else
      {
        CGAL_precondition (equal (max_vertex (cv2), p));
        cv2_to_right = false;        
      }

      // Initialize output flags.
      cv_equal_cv1 = false;
      cv_equal_cv2 = false;
    
      // Take care of the general 4 cases:
      Comparison_result  l_res, r_res;
      Comparison_result  res1, res2;
    
      if (!cv1_to_right && !cv2_to_right)
      {
        // Case 1: Both cv1 and cv2 are defined to the left of p.
        l_res = compare_y_at_x_left (cv1, cv2, p);
      
        if (l_res == LARGER)
        {
          // Case 1(a) : cv1 is above cv2.
          if (!cv_to_right)
          {
            res1 = compare_y_at_x_left (cv1, cv, p);
            res2 = compare_y_at_x_left (cv2, cv, p);
          
            if (res1 == EQUAL)
              cv_equal_cv1 = true;
            if (res2 == EQUAL)
              cv_equal_cv2 = true;
          
            return (res1 == SMALLER || res2 == LARGER);
          }
          return (true);
        }
        else if (l_res == SMALLER)
        {
          // Case 1(b): cv1 is below cv2.
          if (!cv_to_right)
          {
            res1 = compare_y_at_x_left (cv1, cv, p);
            res2 = compare_y_at_x_left (cv2, cv, p);
          
            if (res1 == EQUAL)
              cv_equal_cv1 = true;
            if (res2 == EQUAL)
              cv_equal_cv2 = true;
          
            return (res1 == SMALLER && res2  == LARGER);
          }
          return (false);
        }
        else
        {
          // Overlapping segments.
          if (!cv_to_right)
          {
            res1 = compare_y_at_x_left (cv1, cv, p);
            if (res1 == EQUAL)
            {
              cv_equal_cv1 = true;
              cv_equal_cv2 = true;
              return (false);
            }
            return (true);
          }
          return (true);
        }
      }
      else if (cv1_to_right && cv2_to_right)
      {
        // Case 2: Both cv1 and cv2 are defined to the right of p.
        r_res = compare_y_at_x_right (cv1, cv2, p);

        if (r_res == LARGER)
        {
          // Case 2(a) : cv1 is above cv2.
          if (cv_to_right)
          {
            res1 = compare_y_at_x_right (cv1, cv, p);
            res2 = compare_y_at_x_right (cv2, cv, p);

            if (res1 == EQUAL)
              cv_equal_cv1 = true;
            if (res2 == EQUAL)
              cv_equal_cv2 = true;

            return (res1 == LARGER && res2 == SMALLER);

          }
          return (false);
        }
        else if (r_res == SMALLER)
        {
          // Case 2(b): cv1 is below cv2.
          if (cv_to_right)
          {
            res1 = compare_y_at_x_right (cv1, cv, p);
            res2 = compare_y_at_x_right (cv2, cv, p);

            if (res1 == EQUAL)
              cv_equal_cv1 = true;
            if (res2 == EQUAL)
              cv_equal_cv2 = true;

            return (res1 == LARGER || res2 == SMALLER);
          }
          return (true);
        }
        else
        {
          // Overlapping segments.
          if (cv_to_right)
          {
            res1 = compare_y_at_x_right (cv1, cv, p);
          
            if (res1 == EQUAL)
            {
              cv_equal_cv1 = true;
              cv_equal_cv2 = true;             
              return (false);
            }
            return (true);
          }
          return (true);
        }
      }
      else if (!cv1_to_right && cv2_to_right)
      {
        // Case 3: cv1 is defined to the left of p, and cv2 to its right.
        if (!cv_to_right)
        {
          res1 = compare_y_at_x_left (cv1, cv, p);

          if (res1 == EQUAL)
            cv_equal_cv1 = true;
        
          return (res1 == SMALLER);
        }
        else
        {
          res2 = compare_y_at_x_right (cv2, cv, p);

          if (res2 == EQUAL)
            cv_equal_cv2 = true;

          return (res2 == SMALLER);
        }
      }
      else
      {
        // Case 4: cv1 is defined to the right of p, and cv2 to its left.
        if (cv_to_right)
        {
          res1 = compare_y_at_x_right (cv1, cv, p);

          if (res1 == EQUAL)
            cv_equal_cv1 = true;
        
          return (res1  == LARGER);
        }
        else
        {
          res2 = compare_y_at_x_left (cv2, cv, p);
        
          if (res2 == EQUAL)
            cv_equal_cv2 = true;
 
          return (res2 == LARGER);
        }
      }
    }
  };

  /*! Get an Is_between_cw_2 functor object. */
  Is_between_cw_2 is_between_cw_2_object () const
  {
    return Is_between_cw_2();
  }

  class Compare_cw_around_point_2
  {
  public:
    
    /*!
     * Compare the two interior disjoint x-monotone curves in a clockwise
     * order around their common endpoint.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The common endpoint.
     * \param from_top (true) if we start from 12 o'clock, 
     *                 (false) if we start from 6 o'clock.
     * \pre The point p is an endpoint of both curves.
     * \return SMALLER if we encounter cv1 before cv2;
     *         LARGER if we encounter cv2 before cv1;
     *         EQUAL otherwise.
     */
    Comparison_result operator() (const X_monotone_curve_2& cv1,
                                  const X_monotone_curve_2& cv2,
                                  const Point_2& p,
                                  bool from_top = true) const
    {`aw  q`
      // Find to which side of p (left or right) do cv1 and cv2 lie.
      Base                    tr;
      Construct_min_vertex_2  min_vertex = tr.construct_min_vertex_2_object();
      Construct_max_vertex_2  max_vertex = tr.construct_max_vertex_2_object();
      Equal_2                 equal = tr.equal_2_object();

      bool                    cv1_left, cv2_left;
      
      if (equal (min_vertex (cv1), p)) 
      {
        cv1_left = false;
      }
      else
      {
        CGAL_assertion (equal(max_vertex (cv1), p));
        cv1_left = true;
      }

      if (equal (min_vertex (cv2), p)) 
      {
        cv2_left = false;
      }
      else 
      {
        CGAL_assertion (equal(max_vertex (cv2), p));
        cv2_left = true;
      }

      // Act according to where cv1 and cv2 lie.
      if (cv1_left && cv2_left)
      {
        // Both are defined to the left of p, and we encounter cv1 before
        // cv2 if it is below cv2:
        return (tr.compare_y_at_x_left_2_object() (cv1, cv2, p));
      }
      else if (!cv1_left && !cv2_left)
      {
        // Both are defined to the right of p, and we encounter cv1 before
        // cv2 if it is above cv2. We therefore reverse the order of the
        // curves when we invoke compare_y_at_x_right:
        return (tr.compare_y_at_x_right_2_object() (cv2, cv1, p));
      }
      else if (cv1_left && !cv2_left)
      {
        // If we start from the top, we encounter the right curve (which
        // is cv2) first. If we start from the bottom, we encounter cv1 first.
        if (from_top)
          return (LARGER);
        else
          return (SMALLER);
      }
      else
      {
        // If we start from the top, we encounter the right curve (which
        // is cv1) first. If we start from the bottom, we encounter cv2 first.
        if (from_top)
          return (SMALLER);
        else
          return (LARGER);
      }
    }
  };

  /*! Get a Compare_cw_around_point_2 functor object. */
  Compare_cw_around_point_2 compare_cw_around_point_2_object () const
  {
    return Compare_cw_around_point_2();
  }
  //@}
};

/*! \class
 * A traits-class wrapper that extends the basic traits-class interface.
 */
template <class ArrangementTraits_>
class Arr_traits_wrapper_2 :
  public Arr_traits_basic_wrapper_2<ArrangementTraits_>
{
public:

  // Traits-class geometric types.
  typedef ArrangementTraits_                             Base_traits_2;
  typedef Arr_traits_basic_wrapper_2<ArrangementTraits_> Base;

  typedef typename Base_traits_2::Curve_2                Curve_2;
  typedef typename Base::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Base::Point_2                         Point_2;

  // Tags.
  typedef typename Base::Has_left_category               Has_left_category;
  typedef typename Base::Has_merge_category              Has_merge_category;

  /// \name Construction.
  //@{
  /*! Default constructor. */
  Arr_traits_wrapper_2 () :
    Base()
  {}

  /*! Constructor from a base-traits class. */
  Arr_traits_wrapper_2 (const Base_traits_2& traits) :
    Base (traits)
  {}
  //@}

  // Inherited functors:
  typedef typename Base::Compare_x_2            Compare_x_2;
  typedef typename Base::Compare_xy_2           Compare_xy_2;
  typedef typename Base::Construct_min_vertex_2 Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2 Construct_max_vertex_2;
  typedef typename Base::Is_vertical_2          Is_vertical_2;
  typedef typename Base::Compare_y_at_x_2       Compare_y_at_x_2;
  typedef typename Base::Compare_y_at_x_right_2 Compare_y_at_x_right_2;
  typedef typename Base::Compare_y_at_x_left_2  Compare_y_at_x_left_2;
  typedef typename Base::Equal_2                Equal_2;

  // Note that the basic wrapper does not have to support these functors:
  typedef typename Base_traits_2::Make_x_monotone_2  Make_x_monotone_2;
  typedef typename Base_traits_2::Split_2            Split_2; 
  typedef typename Base_traits_2::Intersect_2        Intersect_2;

  /// \name Overriden functors.
  //@{

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
      // The function is implemented based on the Has_merge category.
      return (_are_mergeable_imp (cv1, cv2, Has_merge_category()));
    }

  private:

    /*!
     * Implementation of the operator() in case the HasMerge tag is true.
     */
    bool _are_mergeable_imp (const X_monotone_curve_2& cv1,
           const X_monotone_curve_2& cv2,
           Tag_true) const
    {
      Base                    tr;
      return (tr.are_mergeable_2_object() (cv1, cv2));      
    }

    /*!
     * Implementation of the operator() in case the HasMerge tag is false.
     */
    bool _are_mergeable_imp (const X_monotone_curve_2& ,
           const X_monotone_curve_2& ,
           Tag_false) const
    {
      // Curve merging is not supported:
      return (false);
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
     *      curve line and share a common endpoint.
     */
    void operator() (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c) const
    {
      // The function is implemented based on the Has_merge category.
      _merge_imp (cv1, cv2, c, Has_merge_category());
    }

  private:

    /*!
     * Implementation of the operator() in case the HasMerge tag is true.
     */
    void _merge_imp (const X_monotone_curve_2& cv1,
                     const X_monotone_curve_2& cv2,
                     X_monotone_curve_2& c,
                     Tag_true) const
    {
      Base                    tr;
      return (tr.merge_2_object() (cv1, cv2, c));      
    }

    /*!
     * Implementation of the operator() in case the HasMerge tag is false.
     */
    void _merge_imp (const X_monotone_curve_2& ,
                     const X_monotone_curve_2& ,
                     X_monotone_curve_2& ,
                     Tag_false) const
    {
      // This function should never be called!
      CGAL_assertion_msg (false,
                          "Merging curves is not supported.");
    }
  };

  /*! Get a Merge_2 functor object. */
  Merge_2 merge_2_object () const
  {
    return Merge_2();
  }
  //@}

};

CGAL_END_NAMESPACE

#endif
