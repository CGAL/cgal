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
// file          : include/CGAL/Planar_map_2/Planar_map_misc_slim.h
// package       : pm (4.08)
// source        :
// revision      :
// revision_date :
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Eyal Flato
//                 Oren Nechushtan <theoren@math.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin)
//
// ======================================================================

#ifndef CGAL_PLANAR_MAP_MISC_H
#define CGAL_PLANAR_MAP_MISC_H

#ifndef CGAL_POLYHEDRON_ITERATOR_3_H
#include <CGAL/Polyhedron_iterator_3.h>
#endif

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
  
  typedef PlanarMapTraits_2          Base;
  typedef typename Base::X_curve     X_curve;
  typedef typename Base::Point_2     Point_2;
  typedef typename Base::Direction_2 Direction_2;
  //typedef  typename PlanarMapTraits_2::Point_2 Point; // for backward compat.
  
  Planar_map_traits_wrap() : Base()
  {
  }

  Planar_map_traits_wrap(const Base& i) : Base(i)
  {
  }

  // Constructing functions
  // ----------------------

  // Access to curve source
  Point_2 curve_source(const X_curve & cv) const 
  { 
    typename Base::Construct_vertex_2 
     construct_vertex = construct_vertex_2_object();
    return construct_vertex(cv, 0);
  }

  // Access to curve target
  Point_2 curve_target(const X_curve & cv) const 
  {
    return construct_vertex_2_object()(cv, 1);
  }

  // Predicate functions
  // -------------------
  // Answers true iff the curve is vertical.
  bool curve_is_vertical(const X_curve & cv) const 
  {
    return is_vertical_2_object()(cv);
  }

  Comparison_result curve_compare_at_x(const X_curve & cv1, 
				       const X_curve & cv2, 
				       const Point_2 & q) const
  {
    return compare_x_at_y_2_objects(q, cv1, cv2);
  }

  bool curve_is_between_cw(const X_curve &cv, 
                           const X_curve &first, 
                           const X_curve &second, 
                           const Point_2 &point) const
  {
    return counterclockwise_in_between_2_object()(cv, first, second, point);
  }
    /*
  {
    Direction_2 d  = construct_direction_2_object()(cv,     point);
    Direction_2 d1 = construct_direction_2_object()(first , point);
    Direction_2 d2 = construct_direction_2_object()(second, point);

    if ( construct_vertex_2_object()(cv, 0)      != point ) 
      d = construct_opposite_direction_2_object()(d);
    if (  construct_vertex_2_object()(first, 0)  != point )
      d1 = construct_opposite_direction_2_object()(d1);
    if (  construct_vertex_2_object()(second, 0) != point )
      d2 = construct_opposite_direction_2_object()(d2);

    return counterclockwise_in_between_2_object()(d, d1, d2);
  }
    */

  // Compares the x value of two points
  Comparison_result compare_x(const Point_2 &p1, const Point_2 &p2) const
  { return compare_x_2_object()(p1, p2); }

  // Compares the y value of two points
  Comparison_result compare_y(const Point_2 &p1, const Point_2 &p2) const
  { return compare_y_2_object()(p1, p2); }

  bool curve_is_same(const X_curve & cv1,const X_curve & cv2) const
  { return equal_2_object()(cv1, cv2); }

  // Intorduce Is_in_x_range_2 / Is_in_x_closed_range_2 ?
  // This can be implemented on the traits_wrap level by using other
  // simpler predicated from the Kernel / traits.
  // Used in the Bounding Box, but there it probably get it from
  // the Pm's Traits_wrap
  bool curve_is_in_x_range(const X_curve & cv, const Point_2 & q) const
  { 
    return !( is_right(q, rightmost(cv.source(), cv.target())) ||
              is_left(q, leftmost(cv.source(), cv.target()))	 );
  }

  /*
  typedef enum
  {
    UNDER_CURVE        = -1,
    CURVE_NOT_IN_RANGE =  0,
    ABOVE_CURVE        =  1,
    ON_CURVE           =  2

  } Curve_point_status;	
  */

  // Returns the curve-point status of the input objects
  typename Base::Curve_point_status 
  curve_get_point_status(const X_curve &cv, const Point_2 & p) const
  {
    if ( ! curve_is_in_x_range(cv, p) )
      return CURVE_NOT_IN_RANGE;
    
    Comparison_result res =
      compare_y_at_x_2_object()(p, cv);
    
    if ( res == LARGER ) return ABOVE_CURVE;
    else if ( res == SMALLER ) return UNDER_CURVE;
    
    return ON_CURVE;
  }
  
  bool point_is_left( const Point_2  & p1, const Point_2  & p2 ) const
  {
    return (compare_x(p1, p2) == SMALLER); 
  }
  
    bool point_is_right( const Point_2  & p1, const  Point_2  & p2 ) const
  {
    return (compare_x(p1, p2) == LARGER); 
  }
  
  bool point_is_same_x( const Point_2  & p1, const Point_2  & p2 ) const
  {
    return (compare_x(p1, p2) == EQUAL); 
  }
  
  bool point_is_lower( const Point_2  & p1, const Point_2  & p2 ) const
  { 
    return (compare_y(p1, p2) == SMALLER); 
  }
  
  bool point_is_higher( const Point_2  & p1, const Point_2  & p2 ) const
  {
    return (compare_y(p1, p2) == LARGER); 
  }
  
  bool point_is_same_y( const Point_2  & p1, const Point_2  & p2 ) const
  {
    return (compare_y(p1, p2) == EQUAL); 
  }
  
  bool point_is_same( const Point_2  & p1, const Point_2  & p2 ) const
  { 
#ifdef PM_MISC_USE_ISSAME
    return is_same(p1, p2);
#else
    return ( (compare_y(p1, p2) == EQUAL) &&
             (compare_x(p1, p2) == EQUAL)   );	
#endif
  }
  
  bool point_is_left_low( const Point_2  & p1,  
                          const Point_2  & p2 ) const
  { 
    Comparison_result k = compare_x(p1, p2);
    if (k == SMALLER)
      return true;
    if ( (k == EQUAL) && (point_is_lower(p1, p2)) )
      return true;
    return false;
  }
  bool point_is_right_top( const Point_2  & p1,  
                          const Point_2  & p2 ) const
  { 
	  return point_is_left_low(p2,p1);
  }
  const Point_2& point_leftmost(const Point_2 &p1, const Point_2 &p2) const
  { return (point_is_left(p1, p2) ? p1 : p2); }
  const Point_2& point_rightmost(const Point_2 &p1, const Point_2 &p2) const
  { return (point_is_right(p1, p2) ? p1 : p2); }
  const Point_2& point_lowest(const Point_2 &p1, const Point_2 &p2) const
  { return (point_is_lower(p1, p2) ? p1 : p2); }
  const Point_2& point_highest(const Point_2 &p1, const Point_2 &p2) const
  { return (point_is_higher(p1, p2) ? p1 : p2); }
  const Point_2& point_leftlow_most(const Point_2 &p1, const Point_2 &p2) const
  { return (point_is_left_low(p1, p2) ? p1 : p2); }
  const Point_2& point_righttop_most(const Point_2 &p1, const Point_2 &p2) const
  { return (point_is_right_top(p1, p2) ? p1 : p2); }
  Point_2 curve_leftmost(const X_curve& cv) const 
  {
    return point_leftmost(curve_source(cv),curve_target(cv));
  }
  Point_2 curve_rightmost(const X_curve& cv) const
  {
    return point_rightmost(curve_source(cv),curve_target(cv));
  }
  Point_2 curve_lowest(const X_curve& cv) const
  {
    return point_lowest(curve_source(cv),curve_target(cv));
  }
  Point_2 curve_highest(const X_curve& cv) const
  {
    return point_highest(curve_source(cv),curve_target(cv));
  }
  Point_2 curve_leftlow_most(const X_curve& cv) const 
  {
    if (!curve_is_vertical(cv)) return curve_leftmost(cv);
    return curve_lowest(cv);
  }
  Point_2 curve_righttop_most(const X_curve& cv) const
  {
    if (!curve_is_vertical(cv)) return curve_rightmost(cv);
    return curve_highest(cv);
  }
  bool curve_merge_condition(const X_curve& whole,
			     const X_curve& part1,
			     const X_curve& part2) const
  {
    return 
      point_is_same(curve_leftlow_most(whole),curve_leftlow_most(part1))&&
      point_is_same(curve_righttop_most(part1),curve_leftlow_most(part2))&&
      point_is_same(curve_righttop_most(whole),curve_righttop_most(part2))||
      point_is_same(curve_leftlow_most(whole),curve_leftlow_most(part2))&&
      point_is_same(curve_righttop_most(part2),curve_leftlow_most(part1))&&
      point_is_same(curve_righttop_most(whole),curve_righttop_most(part1));
  }
  inline bool curve_is_degenerate(const X_curve& cv) const
  {
    return point_is_same(curve_source(cv),curve_target(cv));
  }
public:
  /* precondition:
       cv1,cv2 are adjacent to q
     postcondition:
       returns which of cv1,cv2 is first in clockwise sweep around q
       starting from bottom direction.
  */
  Comparison_result 
  curve_compare_at_x_from_bottom(const X_curve &cv1, 
				 const X_curve &cv2, 
				 const Point_2& q) const 
    {
      if (!curve_is_vertical(cv1))
        if (!curve_is_vertical(cv2))
          if (point_is_same(curve_rightmost(cv1),q)) 
            // cv1 extends leftwards from q
            {
              if (point_is_same(curve_rightmost(cv2),q))
				// cv2 extends leftwards from q
                {
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
				// cv2 extends rightwards from q
                {
                  return curve_compare_at_x_right(cv2,cv1,q);
                }
              else // cv2 extends leftwards from q
                {
                  return LARGER;
                }
            }
        else // cv2 is vertical, cv1 is not vertical 
          {
            if (point_is_same(curve_rightmost(cv1),q) && 
		point_is_same(curve_lowest(cv2),   q))
              return SMALLER;
            else
              return LARGER;
          }
      else // cv1 is vertical
	{
          if (point_is_same(curve_highest(cv1),q))
            if (!curve_is_vertical(cv2) || point_is_same(curve_lowest(cv2),q))
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
  
  Comparison_result 
  curve_compare_at_x_from_top(const X_curve &cv1, 
			      const X_curve &cv2, 
			      const Point_2& q) 
    const 
    {
      if (!curve_is_vertical(cv1))
        if (!curve_is_vertical(cv2))
          if (point_is_same(curve_rightmost(cv1),q)) 
            // cv1 extends leftwards from q
            {
              if (point_is_same(curve_rightmost(cv2),q))
				// cv2 extends leftwards from q
                {
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
				// cv2 extends rightwards from q
                {
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
		point_is_same(curve_highest(cv2), q))
              return SMALLER;
            else
              return LARGER;
          }
      else // cv1 is vertical
	{
          if (point_is_same(curve_lowest(cv1),q))
            if (!curve_is_vertical(cv2) || point_is_same(curve_highest(cv2),q))
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
  bool curve_is_unbounded(const X_curve& cv) const 
  {
    return 
      curve_is_source_unbounded(cv)||
      curve_is_target_unbounded(cv);
  }
};


CGAL_END_NAMESPACE
 

#else   //CGAL_PLANAR_MAP_MISC_H 
#error  Header file Planar_map_misc.h included twice
#endif  //CGAL_PLANAR_MAP_MISC_H 

/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
 *     
 * Planar_map_misc.h - End of File
\*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/
