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

  /*!
   * The function restuns whether two sub-curves can be merged to create
   * the input curve.
   * \param whole The input curve (the merge result).
   * \param part1 The first sub-curve.
   * \param part2 The second sub-curve.
   * \return (true) if whole == part1 + part2.
   */
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
  
  /*!
   * Check whether the curve is degenerate.
   * \param cv The input curve.
   * \return (true) iff the curve source and target are the same.
   */
  inline bool curve_is_degenerate(const X_curve_2 & cv) const
  { 
    return point_is_same(curve_source(cv),curve_target(cv));
  }
    
  /*! 
   * Check if the given query curve is encountered when rotating the first
   * curve in a clockwise direction around a given point until reaching the
   * second curve.
   * \param cv The query curve.
   * \param cv1 The first curve.
   * \param cv2 The second curve.
   * \param p The point around which we rotate cv1.
   * \pre p is an end-point of all three curves.
   * \return (true) if cv is between cv1 and cv2. If cv overlaps cv1 or cv2
   * the result is always (false). If cv1 and cv2 overlap, the result is
   * (true), unless cv1 also overlaps them.
   */
  bool curve_is_between_cw(const X_curve_2& cv, 
                           const X_curve_2& cv1, 
                           const X_curve_2& cv2, 
                           const Point_2& p) const
  {
    // Find the direction of each segment.
    Curve_dir_at_point     dir = _curve_direction_at_point (cv, p);
    Curve_dir_at_point     dir1 = _curve_direction_at_point (cv1, p);
    Curve_dir_at_point     dir2 = _curve_direction_at_point (cv2, p);

    // Special treatment for the cases where cv1 or cv2 are vertical segments:
    if (dir1 == DIR_UP || dir1 == DIR_DOWN)
    {
      if (dir2 == DIR_UP || dir2 == DIR_DOWN)
      {
	// Both cv1 and cv2 are vertical:
	if (dir1 == DIR_UP && dir2 == DIR_DOWN)
	  return (dir == DIR_RIGHT);
	else if (dir1 == DIR_DOWN && dir2 == DIR_UP)
	  return (dir == DIR_LEFT);
	else
	  return (dir != dir1);
      }

      // Only cv1 is vertical:
      if (dir1 == DIR_UP)
      {
	if (dir2 == DIR_LEFT)
	  return (dir == DIR_RIGHT ||
		  dir == DIR_DOWN ||
                  (dir == DIR_LEFT &&
                   curve_compare_at_x_left (cv2, cv, p) == LARGER));
	else
	  return (dir == DIR_RIGHT &&
		  curve_compare_at_x_right (cv2, cv, p) == SMALLER);
      }
      else
      {
	if (dir2 == DIR_LEFT)
	  return (dir == DIR_LEFT &&
		  curve_compare_at_x_left (cv2, cv, p) == LARGER);
	else
	  return (dir == DIR_LEFT ||
		  dir == DIR_UP ||
                  (dir == DIR_RIGHT &&
                   curve_compare_at_x_right (cv2, cv, p) == SMALLER));
      }
    }

    if (dir2 == DIR_UP || dir2 == DIR_DOWN)
    {
      // Only cv2 is vertical:
      if (dir2 == DIR_UP)
      {
	if (dir1 == DIR_LEFT)
	  return (dir == DIR_LEFT &&
		  curve_compare_at_x_left (cv1, cv, p) == SMALLER);
	else
	  return (dir == DIR_LEFT || 
		  dir == DIR_DOWN ||
                  (dir == DIR_RIGHT &&
                   curve_compare_at_x_right (cv1, cv, p) == LARGER));
      }
      else
      {
	if (dir1 == DIR_LEFT)
	  return (dir == DIR_RIGHT ||
		  dir == DIR_UP ||
                  (dir == DIR_LEFT &&
                   curve_compare_at_x_left (cv1, cv, p) == SMALLER));
	else
	  return (dir == DIR_RIGHT &&
		  curve_compare_at_x_right (cv1, cv, p) == LARGER);
      }
    }

    // Take care of the general 4 cases:
    if (dir1 == DIR_LEFT && dir2 == DIR_LEFT)
    {
      // Case 1: Both cv1 and cv2 are defined to the left of p.
      Comparison_result l_res = curve_compare_at_x_left (cv1, cv2, p);
      
      if (l_res == LARGER)
      {
	// Case 1(a) : cv1 is above cv2.
	return (dir != DIR_LEFT ||
		curve_compare_at_x_left (cv1, cv, p) == SMALLER ||
                curve_compare_at_x_left (cv2, cv, p) == LARGER);
      }
      else if (l_res == SMALLER)
      {
	// Case 1(b): cv1 is below cv2.
	return (dir == DIR_LEFT &&
		curve_compare_at_x_left (cv1, cv, p) == SMALLER &&
		curve_compare_at_x_left (cv2, cv, p) == LARGER);
      }
      else
      {
        // Overlapping segments.
        return (dir != DIR_LEFT ||
                curve_compare_at_x_left (cv1, cv, p) != EQUAL);
      }
    }
    else if (dir1 == DIR_RIGHT && dir2 == DIR_RIGHT)
    {
      // Case 2: Both cv1 and cv2 are defined to the right of p.
      Comparison_result r_res = curve_compare_at_x_right (cv1, cv2, p);

      if (r_res == LARGER)
      {
	// Case 2(a) : cv1 is above cv2.
	return (dir == DIR_RIGHT &&
		curve_compare_at_x_right (cv1, cv, p) == LARGER &&
		curve_compare_at_x_right (cv2, cv, p) == SMALLER);
      }
      else if (r_res == SMALLER)
      {
	// Case 2(b): cv1 is below cv2.
	return (dir != DIR_RIGHT ||
		curve_compare_at_x_right (cv1, cv, p) == LARGER ||
                curve_compare_at_x_right (cv2, cv, p) == SMALLER);
      }
      else
      {
        // Overlapping segments.
        return (dir != DIR_RIGHT ||
                curve_compare_at_x_right (cv1, cv, p) != EQUAL);
      }
    }
    else if (dir1 == DIR_LEFT && dir2 == DIR_RIGHT)
    {
      // Case 3: cv1 is defined to the left of p, and cv2 to its right.
      return ((dir == DIR_LEFT &&
	       curve_compare_at_x_left (cv1, cv, p) == SMALLER) ||
	      (dir == DIR_RIGHT &&
	       curve_compare_at_x_right (cv2, cv, p) == SMALLER) ||
	      dir == DIR_UP);
    }
    else
    {
      // Case 4: cv1 is defined to the right of p, and cv2 to its left.
      return ((dir == DIR_RIGHT &&
	       curve_compare_at_x_right (cv1, cv, p) == LARGER) ||
	      (dir == DIR_LEFT &&
	       curve_compare_at_x_left (cv2, cv, p) == LARGER) ||
	      dir == DIR_DOWN);
    }
  }

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
    return (curve_is_source_unbounded(cv)||
	    curve_is_target_unbounded(cv));
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

protected:

  /*!
   * Enum used only be the curve_is_between_cw() function.
   */
  enum Curve_dir_at_point
  {
    DIR_UP,           // Vertical segment, point at 12 o'clock.
    DIR_RIGHT,        // Non-vertical segment going towards the right.
    DIR_DOWN,         // Vertical segment, point at 6 o'clock.
    DIR_LEFT          // Non-vertical segment going towards the left.
  };

  /*!
   * Return the curve direction, with respect to a given refernece point.
   * \param cv The curve.
   * \param p The reference point.
   * \pre p must be an end-point of the segment.
   * \return DIR_UP if cv is a vertical segment pointing at 12 o'clock;
   *         DIR_RIGHT if cv is a non-vertical curve going to the right of p;
   *         DIR_DOWN if cv is a vertical segment pointing at 6 o'clock;
   *         DIR_LEFT if cv is a non-vertical curve going to the left of p;
   */
  Curve_dir_at_point _curve_direction_at_point (const X_curve_2& cv,
						const Point_2& p) const
  {
    // p is one of the end-point. Compare it with the other end-point.
    Comparison_result res;

    if (curve_is_vertical(cv))
    {
      // Special treatment for vertical segments:
      res = compare_xy(p, curve_source(cv));

      if (res == EQUAL)
      {
	res = compare_xy(p, curve_target(cv));
      }
      else
      {
	// Make sure that p is indeed an end-point.
	CGAL_precondition(compare_xy(p, curve_target(cv)) == EQUAL);
      }

      return ((res == SMALLER) ? DIR_UP : DIR_DOWN);
    }

    // In case cv is not vertical:
    res = compare_x(p, curve_source(cv));

    if (res == EQUAL)
    {
      res = compare_x(p, curve_target(cv));
    }
    else
    {
      // Make sure that p is indeed an end-point.
      CGAL_precondition(compare_xy(p, curve_target(cv)) == EQUAL);
    }

    return ((res == SMALLER) ? DIR_RIGHT : DIR_LEFT);
  }
};

CGAL_END_NAMESPACE
 
#endif
