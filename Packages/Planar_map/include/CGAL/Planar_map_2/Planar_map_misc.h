// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Planar_map_2/Planar_map_misc.h
// package       : pm (4.08)
// maintainer    : Eyal Flato           <flato@math.tau.ac.il>
// source        :
// revision      :
// revision_date :
// author(s)     : Iddo Hanniel         <hanniel@math.tau.ac.il>
//                 Eyal Flato           <flato@post.tau.ac.il>
//                 Oren Nechushtan      <theoren@math.tau.ac.il>
//                 Efi Fogel            <efif@post.tau.ac.il>
//                 Ron Wein             <wein@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin)
//
// ======================================================================

#ifndef CGAL_PLANAR_MAP_MISC_H
#define CGAL_PLANAR_MAP_MISC_H

#include <CGAL/config.h>
#include <CGAL/tags.h>

CGAL_BEGIN_NAMESPACE

//--------------------------------------------------------------------------
// Planar_map_traits_wrap - 
//     Geometric Look Up Table. This class extends the user supplied 
// interface to include various less "shallow" operations that are
// impelemented through the interface.
//--------------------------------------------------------------------------

template <class PlanarMapTraits_2>
class Planar_map_traits_wrap : public PlanarMapTraits_2
{
public:
//  typedef  typename I::Info_vertex     Info_vertex;
//  typedef  typename I::Info_edge       Info_edge;
//  typedef  typename I::Info_face       Info_face;
  
  typedef PlanarMapTraits_2                     Base;
  typedef typename Base::X_curve_2              X_curve_2;
  typedef typename Base::Point_2                Point_2;

  typedef typename Base::Has_left_category      Has_left_category;
    
  // Creators:
  // ---------
  Planar_map_traits_wrap() : Base() {}
  Planar_map_traits_wrap(const Base & i) : Base(i) {}

  // Predicates:
  // -----------
  bool point_is_left(const Point_2 & p1, const Point_2 & p2) const
  { return (compare_x(p1, p2) == SMALLER); }
  
  bool point_is_right(const Point_2 & p1, const Point_2 & p2) const
  { return (compare_x(p1, p2) == LARGER); }
    
  bool point_is_same_x(const Point_2 & p1, const Point_2 & p2) const
  { return (compare_x(p1, p2) == EQUAL); }
  
  bool point_is_left_low(const Point_2 & p1, const Point_2 & p2) const
  {
    return (compare_xy(p1, p2) == SMALLER);
  }
    
  bool point_is_right_top(const Point_2 & p1, const Point_2 & p2) const
  { return point_is_left_low(p2,p1); }
    
  const Point_2 & point_leftmost(const Point_2 & p1, const Point_2 & p2) const
  { return (point_is_left(p1, p2) ? p1 : p2); }
    
  const Point_2 & point_rightmost(const Point_2 & p1, const Point_2 & p2) const
  { return (point_is_right(p1, p2) ? p1 : p2); }
    
  const Point_2 & point_leftlow_most (const Point_2 & p1, 
				      const Point_2 & p2) const
  { return (point_is_left_low(p1, p2) ? p1 : p2); }
    
  const Point_2 & point_righttop_most (const Point_2 & p1, 
				       const Point_2 & p2) const
  { return (point_is_right_top(p1, p2) ? p1 : p2); }
    
  Point_2 curve_leftmost(const X_curve_2 & cv) const 
  { return point_leftmost(curve_source(cv),curve_target(cv)); }
    
  Point_2 curve_rightmost(const X_curve_2 & cv) const
  { return point_rightmost(curve_source(cv),curve_target(cv)); }
      
  Point_2 curve_leftlow_most(const X_curve_2 & cv) const 
  {
    if (!curve_is_vertical(cv)) 
      return curve_leftmost(cv);
    else
      return point_leftlow_most(curve_source(cv), curve_target(cv));
  }
    
  Point_2 curve_righttop_most(const X_curve_2 & cv) const
  {
    if (!curve_is_vertical(cv)) 
      return curve_rightmost(cv);
    else
      return point_righttop_most(curve_source(cv), curve_target(cv));
  }
    
  bool curve_merge_condition(const X_curve_2 & whole,
			     const X_curve_2 & part1,
			     const X_curve_2 & part2) const
  {
    // The function simply checks whether it is possible to merge
    // the curves part1 and part2 such that whole is the result.
    if (point_is_same(curve_source(whole), curve_source(part1)))
      if (point_is_same(curve_target(part1), curve_source(part2)))
	if (point_is_same(curve_target(part2), curve_target(whole)))
	  return (true);
	else
	  return (false);
      else
	if (point_is_same(curve_target(part1), curve_target(part2)))
	  if (point_is_same(curve_source(part2), curve_target(whole)))
	    return(true);
	  else
	    return (false);
	else
	  return (false);
    else
      if (point_is_same(curve_source(whole), curve_target(part1)))
	if (point_is_same(curve_source(part1), curve_source(part2)))
	  if (point_is_same(curve_target(part2), curve_target(whole)))
	    return (true);
	  else
	    return (false);
	else
	  if (point_is_same(curve_source(part1), curve_target(part2)))
	    if (point_is_same(curve_source(part2), curve_target(whole)))
	      return (true);
	    else
	      return (false);
	  else
	    return (false);
      else
	if (point_is_same(curve_source(whole), curve_source(part2)) ||
	    point_is_same(curve_source(whole), curve_target(part2)))
	  return (curve_merge_condition (whole, part2, part1));
	else
	  return (false);
  }
    
  inline bool curve_is_degenerate(const X_curve_2 & cv) const
  { return point_is_same(curve_source(cv),curve_target(cv)); }
    
public:
  /*! curve_compare_at_x_from_bottom()
   *
   * \precondition  cv1,cv2 are adjacent to q
   * \postcondition returns which of cv1,cv2 is first in clockwise sweep
   * around q starting from bottom direction.
   */
  Comparison_result 
  curve_compare_at_x_from_bottom(const X_curve_2 & cv1, 
				 const X_curve_2 & cv2, 
				 const Point_2 & q) const
  {
    if (!curve_is_vertical(cv1)) {
      if (!curve_is_vertical(cv2)) {
        if (point_is_same(curve_rightmost(cv1),q))  
        {
          // cv1 extends leftwards from q
          if (point_is_same(curve_rightmost(cv2),q))
          {
            // cv2 extends leftwards from q
            return curve_compare_at_x_left(cv1,cv2,q);
          }
          else // cv2 extends rightwards from q
          {
            return SMALLER;
          }
        }
        else  // cv1 extends rightwards from q
        {
          if (point_is_same(curve_leftmost(cv2),q))
          {
            // cv2 extends rightwards from q
            return curve_compare_at_x_right(cv2,cv1,q);
          }
          else // cv2 extends leftwards from q
          {
            return LARGER;
          }
        }
      } else {
        // cv2 is vertical, cv1 is not vertical 
        if (point_is_same(curve_rightmost(cv1),q) && 
          point_is_same(curve_leftlow_most(cv2),   q))
          return SMALLER;
        else
          return LARGER;
      }
    } else {
      // cv1 is vertical
      if (point_is_same(curve_righttop_most(cv1),q))
        if (!curve_is_vertical(cv2) || 
	    point_is_same(curve_leftlow_most(cv2),q))
          return SMALLER;
        else
          return EQUAL; // both curves extend downwards
      else // cv1 extends from q upwards
        if (point_is_same(curve_righttop_most(cv2),q))
          return LARGER;
        else if (!curve_is_vertical(cv2)) // extends rightwards
          return SMALLER;
        else // cv2 extends upwards
          return EQUAL;
    }
  }

  /*! curve_compare_at_x_from_top()
   */
  Comparison_result 
  curve_compare_at_x_from_top(const X_curve_2 & cv1, 
			      const X_curve_2 & cv2, 
			      const Point_2 & q)
    const 
  {
    if (!curve_is_vertical(cv1))
      if (!curve_is_vertical(cv2))
        if (point_is_same(curve_rightmost(cv1),q)) 
        {
          // cv1 extends leftwards from q
          if (point_is_same(curve_rightmost(cv2),q))
          {
            // cv2 extends leftwards from q
            return curve_compare_at_x_left(cv1,cv2,q);
          }
          else // cv2 extends rightwards from q
          {
            return LARGER;
          }
        }
        else  // cv1 extends rightwards from q
        {
          if (point_is_same(curve_leftmost(cv2),q))
          {
            // cv2 extends rightwards from q
            return curve_compare_at_x_right(cv2,cv1,q);
          }
          else // cv2 extends leftwards from q
          {
            return SMALLER;
          }
        }
      else // cv2 is vertical, cv1 is not vertical 
      {
        if (point_is_same(curve_leftmost(cv1),q) &&
          point_is_same(curve_righttop_most(cv2), q))
          return SMALLER;
        else
          return LARGER;
      }
    else // cv1 is vertical
    {
      if (point_is_same(curve_leftlow_most(cv1),q))
        if (!curve_is_vertical(cv2) || 
	    point_is_same(curve_righttop_most(cv2),q))
          return SMALLER;
        else
          return EQUAL; // both curves extend upwards
      else // cv1 extends from q downwards
        if (point_is_same(curve_leftlow_most(cv2),q))
          return LARGER;
        else if (!curve_is_vertical(cv2)) // extends leftwards
          return SMALLER;
        else // cv2 extends downwards
          return EQUAL;
    }
  }

  /*! curve_is_unbounded()
   */
  bool curve_is_unbounded(const X_curve_2 & cv) const 
  {
    return 
      curve_is_source_unbounded(cv)||
      curve_is_target_unbounded(cv);
  }
    
  /*! curve_compare_at_x_left() is implemented based on the Has_left category
   * If the category indicates that the "left" version is available, it calls
   * the function with same name defined in the base class. Otherwise, it
   * reflects the given point and curves about the origin, and calls the
   * "right" version.
   */
  Comparison_result curve_compare_at_x_left(const X_curve_2 & cv1,
                                            const X_curve_2 & cv2, 
                                            const Point_2 & q) const 
  {
    return curve_compare_at_x_left_imp(cv1, cv2, q, Has_left_category());
  }

    
  Comparison_result curve_compare_at_x_left_imp(const X_curve_2 & cv1,
                                                const X_curve_2 & cv2, 
                                                const Point_2 & q,
                                                Tag_true) const
  {
    return Base::curve_compare_at_x_left(cv1, cv2, q);
  }
    
  Comparison_result curve_compare_at_x_left_imp(const X_curve_2 & cv1,
                                                const X_curve_2 & cv2, 
                                                const Point_2 & q,
                                                Tag_false) const 
  {
    Point_2 rq = point_reflect_in_x_and_y(q);
    X_curve_2 rcv1 = curve_reflect_in_x_and_y(cv1);
    X_curve_2 rcv2 = curve_reflect_in_x_and_y(cv2);
    Comparison_result cr = curve_compare_at_x_right(rcv1, rcv2, rq);
    if (cr == SMALLER) return LARGER;
    if (cr == LARGER) return SMALLER;
    return EQUAL;
  }
};

CGAL_END_NAMESPACE
 
#endif
