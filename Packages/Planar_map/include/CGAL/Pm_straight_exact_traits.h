// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Pm_straight_exact_traits.h
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ============================================================================
#ifndef CGAL_PM_STRAIGHT_EXACT_TRAITS_H
#define CGAL_PM_STRAIGHT_EXACT_TRAITS_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif
#ifndef CGAL_STRAIGHT_2_H
#include <CGAL/Straight_2.h>
#endif

#ifndef CGAL_ISO_RECTANGLE_2_H
#include <CGAL/Iso_rectangle_2.h>
#endif

#ifndef CGAL_INTERSECTION_2_H
#include <CGAL/intersection_2.h>
#endif

#ifndef CGAL_ASSERTIONS_H
#include <CGAL/assertion.h>
#endif

#ifndef CGAL_NUMBER_UTILS_H
#include <CGAL/number_utils.h>
#endif

#ifndef CGAL_PLANAR_MAP_2_ONETUPLE_H
#include <CGAL/Planar_map_2/Onetuple.h>
#endif

#ifndef CGAL_PLANAR_MAP_2_HANDLE_FOR_WITH_ACCESS_H
#include <CGAL/Planar_map_2/Handle_for_with_access.h>
#endif

#ifndef CGAL_PLANAR_MAP_2_NORMALIZE_H
#include <CGAL/Planar_map_2/Normalize.h>
#endif

#ifdef CGAL_PMBB_DEBUG
#include <CGAL/IO/Straight_2_stream.h>
#endif

#define CGAL_PM_STRAIGHT_DEFAULT_NUM 1
#define CGAL_PM_STRAIGHT_DEFAULT_DEN 16
//#define CGAL_PM_STRAIGHT_DEFAULT_DEN 1

CGAL_BEGIN_NAMESPACE

template <class R_>
class Pm_straight_exact_traits
{
public:
  typedef R_ R;
  typedef typename R::FT FT;
  typedef typename R::RT RT;
  typedef Point_2<R> Point;
  typedef Point_2<R> Point_2;
  typedef Vector_2<R> Vector;
  typedef Direction_2<R> Direction;
  typedef Segment_2<R> Segment;
  typedef Ray_2<R> Ray;
  typedef Line_2<R> Line;
  typedef Segment X_bounded_curve;	// [ - curve -> ] 
  typedef Ray X_target_unbounded_curve;	// [ - curve -> ) 
  typedef Ray X_source_unbounded_curve;	// ( - curve -> ]	
  typedef Line X_unbounded_curve;	// ( - curve -> ) 

  // The bounding box is not an Iso rectangle, but a Handle for one
  typedef Iso_rectangle_2<R> Bounding_box;
  typedef Handle_for_with_access<Onetuple<Bounding_box> > 
    Indirect_bounding_box;
  typedef Iso_rectangle_2_Segment_2_pair<R> Bbox_X_bounded_curve;
  typedef Iso_rectangle_2_Ray_2_pair<R> Bbox_X_source_unbounded_curve;
  typedef Iso_rectangle_2_Ray_2_pair<R> Bbox_X_target_unbounded_curve; 
  typedef Iso_rectangle_2_Line_2_pair<R> Bbox_X_unbounded_curve;
  typedef Ray_2_Segment_2_pair<R> Ray_X_bounded_curve;
  typedef Ray_2_Ray_2_pair<R> Ray_X_source_unbounded_curve;
  typedef Ray_2_Ray_2_pair<R> Ray_X_target_unbounded_curve; 
  typedef Ray_2_Line_2_pair<R> Ray_X_unbounded_curve;
  typedef Straight_2_<R> X_curve;
  
  typedef std::vector<Point> Point_container;
  typedef std::vector<X_curve> X_curve_container;
  typedef typename Point_container::iterator Point_iterator;
  typedef typename X_curve_container::iterator X_curve_iterator;
  
  typedef enum 
  {
    CURVE_VERTICAL_UP = 0,
    CURVE_VERTICAL_DOWN = 2,
    CURVE_LEFT = 3,
    CURVE_RIGHT = 1
  } Curve_status;
  
  typedef enum
  {
    UNDER_CURVE = -1,
    ABOVE_CURVE = 1,
    ON_CURVE = 2,
    CURVE_NOT_IN_RANGE = 0
  } Curve_point_status;	

  /*
    typedef enum
    { 
    X_BOUNDED=0,
    X_TARGET_UNBOUNDED=1, 
    X_SOURCE_UNBOUNDED=2, 
    X_UNBOUNDED=8
    } Object_type;
  */	
  typedef enum
  {
    TOP=0,
    RIGHT,
    BOTTOM,
    LEFT
  } Boundary_type;
  
public:
  
  Pm_straight_exact_traits():
    bbox(unbounded_box()),num(CGAL_PM_STRAIGHT_DEFAULT_NUM),
    den(CGAL_PM_STRAIGHT_DEFAULT_DEN)
  {
#ifdef CGAL_PMBB_DEBUG
    std::cout << "\nPm_straight_exact_traits()->";
    debug();
    std::cout << std::flush;
    CGAL_postcondition(is_totally_unbounded());
#endif
  }
  Pm_straight_exact_traits(const RT numer,const RT denum = 1):
    bbox(unbounded_box()),num(numer),den(denum){
    CGAL_precondition(numer>=0 && denum>0);
#ifdef CGAL_PMBB_DEBUG
    std::cout << "\nPm_straight_exact_traits(" << 
      numer << "," << denum << ")->";
    debug();
    std::cout << std::flush;
    CGAL_postcondition(is_totally_unbounded(ref()));
#endif
  }
  Pm_straight_exact_traits(const Bounding_box& b): 
    bbox(b),num(CGAL_PM_STRAIGHT_DEFAULT_NUM),
    den(CGAL_PM_STRAIGHT_DEFAULT_DEN)
  {
#ifdef CGAL_PMBB_DEBUG
    std::cout << 
      "\nPm_straight_exact_traits(";
    debug(b);
    std::cout << ")->";
    debug();
    std::cout << std::flush;
    CGAL_postcondition(bbox==b);
#endif
  }
  Pm_straight_exact_traits(
    const Bounding_box& b, const RT numer,
    const RT denum = 1): 
    bbox(b),num(numer),den(denum){
    CGAL_precondition(numer>=0 && denum>0);
#ifdef CGAL_PMBB_DEBUG
    std::cout << "\nPm_straight_exact_traits(" << b << "," << 
      numer << "," << denum << ")->";
    debug();
    std::cout.flush();
    CGAL_postcondition(b.identical(bbox));
#endif
  }
  // A copy c'tor is required for all the planar map traits.
  Pm_straight_exact_traits(const Pm_straight_exact_traits& tr):
    bbox(tr.bbox),num(tr.num),den(tr.den){
    CGAL_precondition(num>=0 && den>0);
#ifdef CGAL_PMBB_DEBUG
    std::cout << 
      "\nPm_straight_exact_traits(const Pm_straight_exact_traits& tr)->";
    debug();
    std::cout << std::flush;
    CGAL_postcondition(bbox.identical(tr.bbox));
#endif
  }
  inline X_bounded_curve curve_segment(const X_curve & cv) const 
  {
    return curve_segment(cv,ref());
  }
  X_bounded_curve curve_segment(
    const X_curve & cv,
    const Bounding_box& bbox ) const 
  {
    CGAL_precondition_msg(!is_totally_unbounded(bbox),
	"Bounding_box is undefined");

    X_bounded_curve o;
    switch(cv.bound_state())
      {
      case X_curve::NO_UNBOUNDED:
        {
          X_bounded_curve s;
          cv.current(s);
          Bbox_X_bounded_curve bs(&ref(),&s);
          if (!bs.intersection(o))
            {
              CGAL_warning_msg(
                               bs.intersection(o),
                               "\nThe curve must intersect the bounding box");
#ifdef CGAL_PMBB_DEBUG
              debug();
              std::cout << "\ncv " << cv << std::flush;
              std::cout << "\ns " << s << std::flush;
#endif
            }
          break;
        }
      case X_curve::MIN_UNBOUNDED:
        {
          X_source_unbounded_curve r;
          cv.current(r);
          Bbox_X_source_unbounded_curve br(&ref(),&r);
          if (!br.intersection(o))
            {
              CGAL_warning_msg(
                               br.intersection(o),
                               "\nThe curve must intersect the bounding box");
#ifdef CGAL_PMBB_DEBUG
              debug();
              std::cout << "\ncv " << cv << std::flush;
              std::cout << "\nr " << r << std::flush;
#endif
            }
          return o.opposite();
        }
      case X_curve::MAX_UNBOUNDED:
        {
          X_target_unbounded_curve r;
          cv.current(r);
          Bbox_X_target_unbounded_curve br(&ref(),&r);
          if (!br.intersection(o))
            {
#ifdef CGAL_PMBB_DEBUG
              debug();              
              std::cout << "\ncv " << cv << std::flush;
              std::cout << "\nr " << r << std::flush;
#endif
              CGAL_warning_msg(
                               br.intersection(o),
                               "\nThe curve must intersect the bounding box");
            }
          break;
        }
      case X_curve::BOTH_UNBOUNDED:
        {
          X_unbounded_curve l;
          cv.current(l);
          Bbox_X_unbounded_curve bl(&ref(),&l);
          if (!bl.intersection(o))
            {
              CGAL_warning_msg(
                               bl.intersection(o),
                               "\nThe curve must intersect the bounding box");
#ifdef CGAL_PMBB_DEBUG
              debug();
              std::cout << "\ncv " << cv << std::flush;
              std::cout << "\nl " << l << std::flush;
#endif
            }
          break;
        }
      default:
        CGAL_assertion(
                       cv.bound_state()==X_curve::NO_UNBOUNDED||
                       cv.bound_state()==X_curve::MIN_UNBOUNDED||
                       cv.bound_state()==X_curve::MAX_UNBOUNDED||
                       cv.bound_state()==X_curve::BOTH_UNBOUNDED);
      }
    return o;
  }
  
  
  inline Point curve_source(const X_curve & cv,const Bounding_box& bbox) const 
  {
    return curve_segment(cv,bbox).source();
  }
  inline Point curve_source(const X_curve & cv) const 
    // efficient function for calculating the source of the curve
    // Precondition: curve is bounded in the box.
  {
    switch(cv.bound_state())
      {
      case X_curve::NO_UNBOUNDED:
        {
          X_bounded_curve s;
          cv.current(s);
          return s.source();
        }
      case X_curve::MIN_UNBOUNDED:
        {
          X_source_unbounded_curve r;
          cv.current(r);
          return curve_segment(r,ref()).source();
        }
      case X_curve::MAX_UNBOUNDED:
        {
          X_target_unbounded_curve r;
          cv.current(r);
          return r.source();
        }
      case X_curve::BOTH_UNBOUNDED:
        {
          X_unbounded_curve l;
          cv.current(l);
          return curve_segment(cv).source();
        }
#ifdef CGAL_PMBB_DEBUG        
      default:
        CGAL_assertion(bound_state_invariant(cv));
#endif
      }
    return Origin();
  }
  
  inline Point curve_target(const X_curve & cv,const Bounding_box& bbox) const 
  {
    return curve_segment(cv,bbox).target();
  }		
  inline Point curve_target(const X_curve & cv) const 
  {
    switch(cv.bound_state())
      {
      case X_curve::NO_UNBOUNDED:
        {
          X_bounded_curve s;
          cv.current(s);
          return s.target();
        }
      case X_curve::MIN_UNBOUNDED:
        {
          X_source_unbounded_curve r;
          cv.current(r);
          return r.source(); // the target of this ray is actually the source().
        }
      case X_curve::MAX_UNBOUNDED:
        {
          X_target_unbounded_curve r;
          cv.current(r);
          return curve_segment(cv).target();
        }
      case X_curve::BOTH_UNBOUNDED:
        {
          X_unbounded_curve l;
          cv.current(l);
          return curve_segment(cv).target();
        }
#ifdef CGAL_PMBB_DEBUG
      default:
        CGAL_assertion(bound_state_invariant(cv));
#endif
      }
    return Origin();
  }
  
  inline bool curve_is_vertical(const X_curve & cv) const 
  {
    switch (cv.current_state())
      {
      case X_curve::LINE:
        {
          X_unbounded_curve line;
          cv.current(line);
          return line.is_vertical();
        }
      case X_curve::RAY:
        {
          X_target_unbounded_curve ray;
          cv.current(ray);
          return ray.is_vertical();
        }
      case X_curve::SEGMENT:
        {
          X_bounded_curve seg;
          cv.current(seg);
          return seg.is_vertical();
        }
      case X_curve::POINT:
      case X_curve::EMPTY:
        return true;
      }
    CGAL_warning(
                 cv.current_state()==X_curve::LINE||
                 cv.current_state()==X_curve::RAY||
                 cv.current_state()==X_curve::SEGMENT||
                 cv.current_state()==X_curve::POINT||
			  cv.current_state()==X_curve::EMPTY
                 );
    return false;
  }
  
  bool curve_is_in_x_range(const X_curve & cv, const Point & q) const
    { 
      const Point s=curve_source(cv),t=curve_target(cv);
      // precondition: q is inside bbox;cv intersects bbox
      /*
	std::cerr << "\nq=" << q;
	std::cerr << "\ns=" << s;
	std::cerr << "\nt=" << t;
	std::cerr << "\nleftmost(s, t)=" << leftmost(s, t);
	std::cerr << "\nrightmost(s, t)=" << rightmost(s, t);
	std::cerr << "\n!is_right(q, rightmost(s, t))=" << !is_right(q, rightmost(s, t)) ? "true" : "false";
	std::cerr << "\n!is_left(q, leftmost(s, t))=" << !is_left(q, leftmost(s, t)) ? "true" : "false";
	std::cerr.flush();
      */
      return !is_right(q, rightmost(s, t)) && !is_left(q, leftmost(s, t));
    }
  /*
    bool unbounded_curve_is_in_x_range(const X_curve & cv, const Point & q) const
    { 
    switch(cv.get_type())
    {
    case X_BOUNDED:
    return !( is_right(q, rightmost(cv.source(), cv.target())) ||
    is_left(q, leftmost(cv.source(), cv.target()))	 );
    case X_SOURCE_UNBOUNDED:
    return !( is_right(q, rightmost(cv.source(), cv.target())) ||
    is_left(q, leftmost(cv.source(), cv.target()))	 );
    case X_TARGET_UNBOUNDED:
    return !( is_right(q, rightmost(cv.source(), cv.target())) ||
    is_left(q, leftmost(cv.source(), cv.target()))	 );
    case X_UNBOUNDED:
    return !( is_right(q, rightmost(cv.source(), cv.target())) ||
    is_left(q, leftmost(cv.source(), cv.target()))	 );
    
    }
  */
  bool curve_is_in_y_range(const X_curve &cv, const Point & q) const
    { 
      bool r = !( is_lower(q, lowest(curve_source(cv), curve_target(cv))) ||
		  is_higher(q, highest(curve_source(cv), curve_target(cv))) );
      return r;
    }
  
  
  Curve_point_status 
    curve_get_point_status(const X_curve &cv, const Point & p) const
    {
      if (!curve_is_in_x_range(cv, p))
	return CURVE_NOT_IN_RANGE;
      if (!curve_is_vertical(cv))
	{
	  int res = compare_y(p, curve_calc_point(cv, p));
	  if (res == SMALLER) return UNDER_CURVE;
	  if (res == LARGER)	return ABOVE_CURVE;
	  //if (res == EQUAL) 
	  return ON_CURVE;
	}
      else
	{
	  if (is_lower(p,lowest(curve_source(cv),curve_target(cv))))
	    return UNDER_CURVE;
	  if (is_higher(p,highest(curve_source(cv),curve_target(cv))))
	    return ABOVE_CURVE;
	  // if (curve_is_in_y_range(cv,p))
	  return ON_CURVE;
	}
    }
  /*
    Curve_point_status 
    unbounded_curve_get_point_status(const X_curve &cv, const Point & p) const
    {
    if(cv.is_vertical())
    {
    if (is_same_x(cv.point(),p)) return 
    
    if (!unbounded_curve_is_in_x_range(cv, p))
    return CURVE_NOT_IN_RANGE;
    if (!curve_is_vertical(cv))
    {
    int res = compare_y(p, curve_calc_point(cv, p));
    if (res == SMALLER) return UNDER_CURVE;
    if (res == LARGER)	return ABOVE_CURVE;
    //if (res == EQUAL) 
    return ON_CURVE;
    }
    else
    {
    if (is_lower(p,lowest(curve_source(cv),curve_target(cv))))
    return UNDER_CURVE;
    if (is_higher(p,highest(curve_source(cv),curve_target(cv))))
    return ABOVE_CURVE;
    // if (curve_is_in_y_range(cv,p))
    return ON_CURVE;
    }
    }
  */	  
  Comparison_result 
    curve_compare_at_x(const X_curve &cv1, const X_curve &cv2, const Point &q) 
    const 
    {
      //CGAL_assertion (curve_is_in_x_range(cv1, q));
      //CGAL_assertion (curve_is_in_x_range(cv2, q));
      if ((!curve_is_in_x_range(cv1, q)) || (!curve_is_in_x_range(cv2, q)))
	return EQUAL;
      
      Point p1 = curve_calc_point(cv1, q);
      Point p2 = curve_calc_point(cv2, q);
      
      if (curve_is_vertical(cv1))
	{
	  if (curve_is_vertical(cv2))
	    {
	      // both cv1 and cv2 are vertical
	      if ( is_lower(curve_target(cv1), curve_source(cv2)) )
		return SMALLER;
	      if ( is_higher(curve_source(cv1), curve_target(cv2)) )
		return LARGER;
	      return SMALLER;
	    }
	  // cv1 is vertical and cv2 not
	  if ( is_lower(curve_target(cv1), p2) )
	    return SMALLER;
	  if ( is_higher(curve_source(cv1), p2) )
	    return LARGER;
	  return EQUAL;
	}
      
      if (curve_is_vertical(cv2))
	{
	  // cv2 is vertical and cv1- not
	  /*        bug fix (Oren)
		    if (is_lower(curve_target(cv2), p1) )
		    return LARGER;
		    if ( is_higher(curve_source(cv2), p1) )
		    return SMALLER;
		    
		    if ( is_higher(curve_source(cv2), p1) ) // bug fix (Oren)
		    The answer should be independent of the curve's orientation !!
		    
		    p1 x--x               p1 x--x
		    |                     /\
		    |           versus    |
		    \/cv2                 |cv2
		    x                     x
		    
		    p                     p
		    
	  */
	  if (is_lower(lowest(curve_source(cv2),curve_target(cv2)), p1) )
	    return LARGER;
	  if ( is_higher(highest(curve_source(cv2),curve_target(cv2)), p1) )
	    return SMALLER;
	  return EQUAL;
	}
      
      // both are not vertical
      if (is_higher(p1, p2)) return LARGER;
      if (is_lower(p1, p2)) return SMALLER;
      return EQUAL;
    }
  
  
  Comparison_result 
    curve_compare_at_x_left(const X_curve &cv1, const X_curve &cv2, 
			    const Point &q) const 
    {
      // cases  in which the function isn't defined
      //CGAL_assertion(!curve_is_vertical(cv1));
      //CGAL_assertion(!curve_is_vertical(cv2));
		  //CGAL_assertion(is_left(leftmost(curve_source(cv1), curve_target(cv1)), q));
		  //CGAL_assertion(is_left(leftmost(curve_source(cv2), curve_target(cv2)), q));
		  
		  if (curve_is_vertical(cv1) || (curve_is_vertical(cv2))) return EQUAL;
		  if (!is_left(leftmost(curve_source(cv1), curve_target(cv1)), q)) return EQUAL;
		  if (!is_left(leftmost(curve_source(cv2), curve_target(cv2)), q)) return EQUAL;
		  
		  Comparison_result r = curve_compare_at_x(cv1, cv2, q);
		  
		  if ( r != EQUAL)
			  return r;     // since the curve is continous 
		  
		  // <cv2> and <cv1> meet at a point with the same x-coordinate as q
		  // compare their derivatives
		  return compare_derivative(cv2,cv1);
	  }
	  
  Comparison_result 
    curve_compare_at_x_right(const X_curve &cv1, const X_curve &cv2, const Point & q) const 
    {
      // cases  in which the function isn't defined
      //CGAL_assertion(!curve_is_vertical(cv1));
      //CGAL_assertion(!curve_is_vertical(cv2));
      //CGAL_assertion(is_right(rightmost(curve_source(cv1), curve_target(cv1)), q));
      //CGAL_assertion(is_right(rightmost(curve_source(cv2), curve_target(cv2)), q));
      
      if (curve_is_vertical(cv1) || (curve_is_vertical(cv2))) return EQUAL;
      if (!is_right(rightmost(curve_source(cv1), curve_target(cv1)), q)) return EQUAL;
      if (!is_right(rightmost(curve_source(cv2), curve_target(cv2)), q)) return EQUAL;
      
      Comparison_result r = curve_compare_at_x(cv1, cv2, q);
      
      if ( r != EQUAL)
	return r;     // since the curve is continous (?)
      
      // <cv1> and <cv2> meet at a point with the same x-coordinate as q
      // compare their derivatives
      return compare_derivative(cv1,cv2);
    }
  
  
  const X_curve curve_flip(const X_curve &cv) const
    {
      switch(cv.bound_state())
	{
	case X_curve::NO_UNBOUNDED:
	  {
	    X_bounded_curve seg;
	    cv.current(seg);
	    return seg.opposite();
	  }
	  // same curve with opposite orientation.
	case X_curve::MIN_UNBOUNDED:
	  {
	    X_target_unbounded_curve ray;
	    cv.current(ray);
	    return X_curve(ray,true);
	  }			  
	  // same curve with opposite orientation.
	case X_curve::MAX_UNBOUNDED:
	  {
	    X_target_unbounded_curve ray;
	    cv.current(ray);
	    return X_curve(ray,false);
	  }			  
	case X_curve::BOTH_UNBOUNDED:
	  {
	    X_unbounded_curve line;
	    cv.current(line);
	    return line.opposite();
	  }		  
	default:
	  CGAL_assertion(
			 cv.bound_state()==X_curve::NO_UNBOUNDED||
			 cv.bound_state()==X_curve::MIN_UNBOUNDED||
			 cv.bound_state()==X_curve::MAX_UNBOUNDED||
			 cv.bound_state()==X_curve::BOTH_UNBOUNDED);
	}
      return X_curve();
    }
/*
  const Direction direction(const X_curve &cv) const
  {
    switch(cv.current_state())
      {
      case X_BOUNDED:
        return ((X_bounded_curve&)cv).direction();
      case X_TARGET_UNBOUNDED:
      case X_SOURCE_UNBOUNDED:
        return ((X_target_unbounded_curve&)cv).direction();
      case X_UNBOUNDED:
        return ((X_unbounded_curve&)cv).direction();
      default:
        CGAL_assertion(cv.current_state()==X_BOUNDED||cv.current_state()==X_TARGET_UNBOUNDED||
                       cv.current_state()==X_SOURCE_UNBOUNDED||cv.current_state()==X_UNBOUNDED);
      }
    return Direction();
  }
*/	  
  Curve_status curve_get_status(const X_curve &cv) const
    {
      if (curve_is_vertical(cv)) 
	{
	  if ( is_higher(curve_target(cv), curve_source(cv)) )
	    return CURVE_VERTICAL_UP;
	  else
	    return CURVE_VERTICAL_DOWN;
	}
      else
	{
	  if ( is_right(curve_target(cv), curve_source(cv)) )
	    return CURVE_RIGHT;
	  else
	    return CURVE_LEFT;
	}
    }
  
  bool curve_is_between_cw(const X_curve &cv, 
			   const X_curve &first, 
			   const X_curve &second, 
			   const Point &cp)	const

    // TRUE if cv is between first and second in cw direction
    // precondition: cv, first and second have a common endpoint
    // precondition: first, second, cv are pairwise interior disjoint
    {
#ifdef CGAL_PMBB_DEBUG
      std::cerr << "\nbool curve_is_between_cw(" << cv << "," << first << "," 
                << second << "," << cp << ")" << std::endl;
#endif
      // CGAL_assertion(is_intersection_simple(first, second);
      // CGAL_assertion(is_intersection_simple(first, *this);
      // CGAL_assertion(is_intersection_simple(*this, second);
      
      Curve_status cv0_status, cv1_status, cvx_status;
      int cv0_cv1, cv0_cvx, cv1_cvx;
      cv0_cv1 = cv0_cvx = cv1_cvx = -1;
      
      X_curve cv0 = first;
      X_curve cv1 = second;
      X_curve cvx = cv;
      
      if ( !is_same(curve_source(cv0),cp) ) cv0 = curve_flip(cv0);
      if ( !is_same(curve_source(cv1),cp) ) cv1 = curve_flip(cv1);
      if ( !is_same(curve_source(cvx),cp) ) cvx = curve_flip(cvx);
      
      cv0_status = curve_get_status(cv0);
      cv1_status = curve_get_status(cv1);
      cvx_status = curve_get_status(cvx);
      
      //	the circle:				    0
      //						 ** | **
      //						*	*
      //					     3 *	 * 1
      //						*	*
      //						 ** | **
      //						    2
      
      if (cv0_status == cv1_status)
	{
	  if (cv0_status == CURVE_RIGHT)
	    cv0_cv1 = curve_compare_at_x_right(cv0, cv1, cp);
	  else if (cv0_status == CURVE_LEFT)
	    cv0_cv1 = curve_compare_at_x_left(cv0, cv1, cp);
	}
      if (cv0_status == cvx_status)
	{
	  if (cv0_status == CURVE_RIGHT)
	    cv0_cvx = curve_compare_at_x_right(cv0, cvx, cp);
	  else if (cv0_status == CURVE_LEFT)
	    cv0_cvx = curve_compare_at_x_left(cv0, cvx, cp);
	}
      if (cv1_status == cvx_status)
	{
	  if (cv1_status == CURVE_RIGHT)
	    cv1_cvx = curve_compare_at_x_right(cv1, cvx, cp);
	  if (cv1_status == CURVE_LEFT)
	    cv1_cvx = curve_compare_at_x_left(cv1, cvx, cp);
	}
      
      if (cv0_status == cv1_status)
	{
	  if (cv0_status == CURVE_LEFT)
	    {
	      if ( ((cv0_cv1==1) && (cvx_status==cv0_status) && 
		    ((cv0_cvx==-1) || (cv1_cvx==1))) ||
		   ((cv0_cv1==1) && (cvx_status!=cv0_status)) ||
		   ((cv0_cv1==-1) && (cvx_status==cv0_status) && 
		    ((cv0_cvx==-1) && (cv1_cvx==1))) )
		return true;
	    }
	  if (cv0_status == CURVE_RIGHT)
	    {
	      if ( ((cv0_cv1==1) && (cvx_status==cv0_status) && 
		    ((cv0_cvx==1) && (cv1_cvx==-1))) ||
		   ((cv0_cv1==-1) && (cvx_status!=cv0_status)) ||
		   ((cv0_cv1==-1) && (cvx_status==cv0_status) && 
		    ((cv0_cvx==1) || (cv1_cvx==-1))) )
		return true;
	    }
	  return false;
	}
      // else do the following
      
      if (cv0_status == cvx_status)
	{
	  if ( ((cv0_status == CURVE_LEFT) && (cv0_cvx==-1)) ||
	       ((cv0_status == CURVE_RIGHT) && (cv0_cvx==1)) )
	    return true;
	  
	  //Addition by iddo for enabeling addition of null segments - testing
	  if ( (cv0_status==CURVE_VERTICAL_DOWN)&&
	       ((curve_source(cv0)==curve_target(cv0))||
		(curve_source(cvx)==curve_target(cvx))) )
	    return true; //a null segment (=point) 
	  
	  return false;
	}
      
      if (cv1_status == cvx_status)
	{
	  if ( ((cv1_status == CURVE_LEFT) && (cv1_cvx==1)) ||
	       ((cv1_status == CURVE_RIGHT) && (cv1_cvx==-1)) )
	    return true;
	  
	  //Addition by iddo for enabeling addition of null segments - testing
	  if ( (cv1_status==CURVE_VERTICAL_DOWN)&&
	       ((curve_source(cv1)==curve_target(cv1))
		||(curve_source(cvx)==curve_target(cvx))) )
	    return true; //a null segment (=point)  
	  
	  return false;
	}
      
      // cv1 and cv0 are on diffrent part of the circle - it is easy
      if ( ((cv1_status - cv0_status + 4)%4) < 
	   ((cvx_status - cv0_status + 4)%4) )
	return false;
      else
	// if there is an equality or inequality to the other side
	// everything is ok
	return true;
    }
  
  Comparison_result compare_x(const Point &p1, const Point &p2) const
    { return compare_value(p1.x(), p2.x()); }
  Comparison_result compare_y(const Point &p1, const Point &p2) const
    { return compare_value(p1.y(), p2.y()); }
 public:
  Point point_to_left(const Point& p) const {
    return p+Vector(-num,0,denum);
  }
  Point point_to_right(const Point& p) const {
    return p+Vector(num,0,denum);
  }
  bool curve_is_same(const X_curve &cv1, const X_curve &cv2) const
    {
      return cv1==cv2;
    }
  bool bounding_box_is_same(const Bounding_box& b) const{
/* returns true if trait's bounding box and input's are the same */
    if (is_totally_unbounded()) return is_totally_unbounded(b);
    if (is_totally_unbounded(b)) return false;
// NOTE: operator== for CGAL homo. points returns true if one of the
// points is homogeneous with last entry 0, that is, infinite.
    return ref()==b;
  }
  bool is_point_on_curve(const X_curve &cv, const Point& p) const //check
    {
      if (!curve_is_in_x_range(cv, p))
	return false;
      if (curve_is_vertical(cv))
	{
	  if (curve_is_in_y_range(cv,p))
	    return true;
	  else
	    return false;
	}
      int res = compare_y(p, curve_calc_point(cv, p));
      if (res == EQUAL)
	return true;
      return false;
    }
 private:
  bool is_left(const Point &p1, const Point &p2) const 
    { return (compare_x(p1, p2) == SMALLER); }
	bool is_right(const Point &p1, const Point &p2) const 
	{ return (compare_x(p1, p2) == LARGER); }
	bool is_same_x(const Point &p1, const Point &p2) const 
	{ return (compare_x(p1, p2) == EQUAL); }
	bool is_lower(const Point &p1, const Point &p2) const 
	{ return (compare_y(p1, p2) == SMALLER); }
	bool is_higher(const Point &p1, const Point &p2) const 
	{ return (compare_y(p1, p2) == LARGER); }
	bool is_same_y(const Point &p1, const Point &p2) const 
	{ return (compare_y(p1, p2) == EQUAL); }
	bool is_same(const Point &p1, const Point &p2) const
	{
		return (compare_x(p1, p2) == EQUAL) &&
			(compare_y(p1, p2) == EQUAL);
	}
	const Point& leftmost(const Point &p1, const Point &p2) const
	{ return (is_left(p1, p2) ? p1 : p2); }
	
	const Point& rightmost(const Point &p1, const Point &p2) const
	{ return (is_right(p1, p2) ? p1 : p2); }
	
	const Point& lowest(const Point &p1, const Point &p2) const
	{ return (is_lower(p1, p2) ? p1 : p2); }
	
	const Point& highest(const Point &p1, const Point &p2) const
	{ return (is_higher(p1, p2) ? p1 : p2); }
	
	bool is_left(const X_curve&cv, const Point &p) const 
	{ return (compare_x(rightmost(cv), p) == SMALLER); }
	bool is_right(const X_curve &cv, const Point &p) const 
	{ return (compare_x(leftmost(cv), p) == LARGER); }
	const Point leftmost(const X_curve& cv) const 
	{
		return leftmost(curve_source(cv),curve_target(cv));
	}
	const Point rightmost(const X_curve& cv) const
	{
		return rightmost(curve_source(cv),curve_target(cv));
	}
	
	const Point lowest(const X_curve& cv) const
	{
		return lowest(curve_source(cv),curve_target(cv));
	}
	const Point highest(const X_curve& cv) const
	{
		return highest(curve_source(cv),curve_target(cv));
	}
	/*	
	const Point lexleftmost(const X_curve& cv) const
	{ 	  
	if (!curve_is_vertical(cv)) return leftmost(cv);
	return lowest(cv);
	}
	const Point lexrightmost(const X_curve& cv) const
	{ 
	if (!curve_is_vertical(cv)) return rightmost(cv);
	return highest(cv);
	}	
	*/
	const Point& lexleftmost(const Point& p,const Point& q) const
	{ 	  
		if (!is_same_x(p,q)) return leftmost(p,q);
		return lowest(p,q);
	}
	const Point& lexrightmost(const Point& p,const Point& q) const
	{ 
		if (!is_same_x(p,q)) return rightmost(p,q);
		return highest(p,q);
	}	
	
public:
	Point curve_calc_point_old(const X_curve &cv, const Point & q) const
		// 	Used to draw an arrow representation of the vertical ray shoot.
	{
		// CGAL_assertion (!curve_is_in_s_range(cv, q));
		if ( !curve_is_in_x_range(cv, q) )
			return curve_source(cv);
		
		if (curve_is_vertical(cv))
			return curve_source(cv);
		
		//return Point(q.x(), curve_source(cv).y() + 
		//             (curve_target(cv).y() - curve_source(cv).y()) / 
		//             (curve_target(cv).x() - curve_source(cv).x()) * 
		//             (q.x() - curve_source(cv).x()) );
		
		const Point & a = curve_source(cv);
		const Point & b = curve_target(cv);
		return Point ((b.hx() * a.hw() - a.hx() * b.hw()) * q.hx() * a.hw(),
			(b.hx() * a.hw() - a.hx() * b.hw()) * q.hw() * a.hy() + 
			(b.hy() * a.hw() - a.hy() * b.hw()) * 
			(q.hx() * a.hw() - a.hx() * q.hw()),  
			(b.hx() * a.hw() - a.hx() * b.hw()) * q.hw() * a.hw());
	}
  Point curve_calc_point(const X_curve &cv, const Point & q) const
    // 	Used to draw an arrow representation of the vertical 
    // ray shoot.
  {
    const Point s=curve_source(cv),t=curve_target(cv);
    if (is_right(t,s))
      {
        if (is_right(q,t)) return t;
        if (is_left(q,s)) return s;
      }
    else if (is_left(t,s))
      {
        if (is_right(q,s)) return s;
        if (is_left(q,t)) return t;
      }
    else // is_same_x(t,s)
      return s;
    /* returns any point on the segment. */
    
    //return Point(q.x(), curve_source(cv).y() + 
    //             (curve_target(cv).y() - curve_source(cv).y()) / 
    //             (curve_target(cv).x() - curve_source(cv).x()) * 
    //             (q.x() - curve_source(cv).x()) );
    Point res1((t.hx() * s.hw() - s.hx() * t.hw()) * q.hx() * s.hw(),
               //			(t.hx() * s.hw() - s.hx() * t.hw()) * q.hw() * s.hy() + 
               //			(t.hy() * s.hw() - s.hy() * t.hw()) * 
               //			(q.hx() * s.hw() - s.hx() * q.hw()),  
               (t.hx() * q.hw() - q.hx() * t.hw()) * s.hw() * s.hy() + 
               (q.hx() * s.hw() - s.hx() * q.hw()) * t.hy() * s.hw(),  
               (t.hx() * s.hw() - s.hx() * t.hw()) * q.hw() * s.hw());
    
    Point res2((t.hx() * s.hw() - s.hx() * t.hw()) * q.hx() * s.hw(),
               (t.hx() * s.hw() - s.hx() * t.hw()) * q.hw() * s.hy() + 
               (t.hy() * s.hw() - s.hy() * t.hw()) * 
               (q.hx() * s.hw() - s.hx() * q.hw()),  
               (t.hx() * s.hw() - s.hx() * t.hw()) * q.hw() * s.hw());
    
#ifdef CGAL_PMBB_DEBUG
    switch (cv.current_state())
      {
      case X_curve::EMPTY:
      case X_curve::POINT:
        break;
      case X_curve::LINE:
        {
          X_unbounded_curve line;
          cv.current(line);
          CGAL_warning(line.has_on_boundary(res1));
          CGAL_warning(line.has_on_boundary(res2));
          break;
        }
      case X_curve::RAY:
        {
          X_target_unbounded_curve ray;
          cv.current(ray);
          CGAL_warning(ray.has_on(res1));
          CGAL_warning(ray.has_on(res2));
          CGAL_warning(ray.collinear_has_on(res1));
          CGAL_warning(ray.collinear_has_on(res2));
          break;
        }
      case X_curve::SEGMENT:
        {
          X_bounded_curve seg;
          cv.current(seg);
          //CGAL_warning(seg.has_on(res1)); // very sensitive function.
          //CGAL_warning(seg.has_on(res2)); // e.g., when using doubles.
          CGAL_warning(seg.collinear_has_on(res1));
          CGAL_warning(seg.collinear_has_on(res2));					
          break;
        }
      }
#endif
    return res1;
  }
  
private:	
  Comparison_result compare_derivative(const X_curve &cv1,const X_curve& cv2) const
  {
    FT dx1=curve_target(cv1).x() - curve_source(cv1).x(),
      dx2=curve_target(cv2).x() - curve_source(cv2).x(),
      dy1=curve_target(cv1).y() - curve_source(cv1).y(),
      dy2=curve_target(cv2).y() - curve_source(cv2).y(),zero(0);
    return static_cast<Comparison_result>(compare_value(dx1,zero)*compare_value(dx2,zero)*compare_value(dy1*dx2,dy2*dx1));
  }
  /*  
  typename R::FT curve_b_const(const X_curve &cv)const
  {
    CGAL_assertion (!curve_is_vertical(cv));
    return ((curve_target(cv).x() * curve_source(cv).y() - 
             curve_target(cv).y()*curve_source(cv).x())     / 
            (curve_target(cv).x() - curve_source(cv).x()));
  }
  */
  Comparison_result compare_value(const typename R::FT &v1, 
                                  const typename R::FT &v2) const
    {
      typename R::FT d = v1 - v2;
      typename R::FT z(0);
      if (d == z)
	return EQUAL;
      if (z < d)
	return LARGER;
      else
	return SMALLER;
    }
  
 public:
  inline const Bounding_box& get_bounding_box() const {return ref();}
  /* Returns the bbox's point boundary starting from left bottom 
     counter clockwise. */
protected:
  inline const Point get_point_boundary(const Boundary_type& i,const Bounding_box& b) const {
    CGAL_precondition_msg(!is_totally_unbounded(b),"Bounding_box is undefined");
		return b[i];
  }
  inline const Point get_point_boundary(const Boundary_type& i) const {
    return get_point_boundary(i,ref());
  }
public:
  inline void get_point_boundary(Point_container& c,const Bounding_box& b) {
    for (int i=0;i<4;i++) c[i]=get_point_boundary(Boundary_type(i),b);
	}
  inline void get_point_boundary(Point_container& c) {
    for (int i=0;i<4;i++) c[i]=get_point_boundary(Boundary_type(i),ref());
  }
protected:
  inline const X_curve get_x_curve_boundary(const Boundary_type& i,
                                            const Bounding_box& b) const {
    /* returns the bounding box x_curve boundary starting from 
       bottom counter clockwise. */
    CGAL_precondition_msg(!is_totally_unbounded(b),"Bounding_box is undefined");
		CGAL_warning_msg(i>=0 && i<4,"\nVertex should be within the range {0,1,2,3}");
		return Segment(b[i],b[(i+1)%4]);
  }
  inline const X_curve get_x_curve_boundary(const Boundary_type& i) const {
    return get_x_curve_boundary(i,ref());
  }
public:
  inline void get_x_curve_boundary(X_curve_container& c,const Bounding_box& b)
    const {
		for (int i=0;i<4;i++)
                  c.push_back(X_curve(Segment(b[i],b[(i+1)%4])));
  }
  inline void get_x_curve_boundary(X_curve_container& c) const {
    get_x_curve_boundary(c,ref());
  }
  bool curve_is_source_unbounded(const X_curve& cv) const{
    switch (cv.bound_state())
      {
      case X_curve::NO_UNBOUNDED:
      case X_curve::MAX_UNBOUNDED:
      case X_curve::LINE_EMPTY:
        return false;
      case X_curve::MIN_UNBOUNDED:
      case X_curve::BOTH_UNBOUNDED:
        return true;
      }
#ifdef CGAL_PMBB_DEBUG
    debug_invariant(bound_state_invariant(cv));
#endif
    return false;
  }

  bool curve_is_target_unbounded(const X_curve& cv) const{
    // We require this additional function for efficiency sake.
    switch (cv.bound_state())
      {
      case X_curve::NO_UNBOUNDED:
      case X_curve::MIN_UNBOUNDED:
      case X_curve::LINE_EMPTY:
	return false;
      case X_curve::MAX_UNBOUNDED:
      case X_curve::BOTH_UNBOUNDED:
	return true;
      }
#ifdef CGAL_PMBB_DEBUG
    debug_invariant(bound_state_invariant(cv));
#endif
    return false;
  }
  
  Bounding_box get_bounding_box(const Point& p) const
    /* Returns a bounding box that encloses the point */
    {
      Bounding_box bb=Bounding_box(p+Vector(-num,-num,den),
				   p+Vector(num,num,den));
      normalize_coordinates(bb);
      /*
#ifdef CGAL_PMBB_DEBUG
      debug_invariant(is_bounded_invariant(p,bb));
#endif
      */
      return bb;
    }
  
  Bounding_box get_bounding_box(const X_curve& cv) const {
    /* Return a bounding box that encloses the curve */
    const Bounding_box& bbox=ref();
      switch(cv.current_state())
	{
	case X_curve::POINT:
	  return bbox;
	  /*		
	    currently inserting degenerate curves are not supported  
	    {
	    Point p;
	    cv.current(p);
	    Bounding_box bb=Bounding_box(p+Vector(-num,-num,den),
	                                 p+Vector(num,num,den));
	    normalize_coordinates(bb);
	    return bb;
	    }
	  */
	case X_curve::SEGMENT:
	  {
	    X_bounded_curve s;
	    cv.current(s);
	    const Bounding_box curr(s.source(),s.target());
	    Bounding_box bb = Bounding_box(
		       curr[0]+Vector(-num,-num,den),	//bbox[0] - left bottom
	               curr[2]+Vector(num,num,den));	//bbox[2] - right top
	    normalize_coordinates(bb);
            /*
#ifdef CGAL_PMBB_DEBUG
            debug_invariant(is_bounded_invariant(s,bb));
#endif
            */
	    return bb;
	  }
	case X_curve::RAY:
	  {
	    X_target_unbounded_curve r;
	    cv.current(r);
	    const Point& p=r.source(); 
				// Ray is always X_target_unbounded_curve
	    Bounding_box bb=Bounding_box(
	      p+Vector(-num,-num,den),
	      p+Vector(num,num,den));
	    normalize_coordinates(bb);
            /*
#ifdef CGAL_PMBB_DEBUG
      debug_invariant(is_bounded_invariant(r,bb));
#endif
            */
	    return bb;
	  }
	case X_curve::LINE:
	  {
	    X_unbounded_curve l;
	    cv.current(l);
	    const Point& p=l.point();
	    Bounding_box bb=Bounding_box(
	       p+Vector(-num,-num,den),
	       p+Vector(num,num,den));
	    normalize_coordinates(bb);
	    return bb;
	  }
	case X_curve::EMPTY:
	  return bbox;
	}	
      CGAL_assertion_msg(	
	cv.current_state()==X_curve::POINT||
	cv.current_state()==X_curve::SEGMENT||
	cv.current_state()==X_curve::RAY||
	cv.current_state()==X_curve::LINE||
	cv.current_state()==X_curve::EMPTY,
	"Wrong curve type in\nconst Bounding_box get_bounding_box(const X_curve& cv) const");
      return bbox;
    }

  const Bounding_box get_bounding_box(const Point& p,
				      const Bounding_box& b) const
    {
      CGAL_precondition_msg(!is_totally_unbounded(b),
			    "Bounding_box is undefined");
      if (is_point_bounded(p,b)) return b;
      Bounding_box out = Bounding_box(
	Bounding_box(b[0],p+Vector(-num,-num,den))[0],
	Bounding_box(b[2],p+Vector(num,num,den))[2]);
      normalize_coordinates(out);
      CGAL_postcondition_msg(!is_totally_unbounded(out),
			     "Bounding_box is undefined");
#ifdef CGAL_PMBB_DEBUG
      debug_invariant(is_point_bounded(p,out));
      debug_invariant(is_bounding_box_bounded(b,out));
#endif
      return out;
    }

  const Bounding_box get_bounding_box(const X_curve& cv,
                                      const Bounding_box& b) const
/* Returns a bounding box that encloses both the curve and the box */
  {
    CGAL_precondition_msg(!is_totally_unbounded(b),
                          "Bounding_box is undefined");
    if (is_curve_bounded(cv,b)) return b;
    Bounding_box out = ref();
    switch(cv.current_state())
      {
      case X_curve::EMPTY:
      case X_curve::POINT: // currently not supported
        return out;
      case X_curve::SEGMENT:
        {
          if (is_curve_bounded(cv,b)) return b;
          X_bounded_curve s;
          cv.current(s);
          const Bounding_box curr(s.source(),s.target());
          out=Bounding_box(
                   Bounding_box(b[0],curr[0]+Vector(-num,-num,den))[0],
                   //bbox[0] - left bottom
                   Bounding_box(b[2],curr[2]+Vector(num,num,den))[2]);	
                   //bbox[2] - right top
	  normalize_coordinates(out);
          break;
        }
      case X_curve::RAY:
        {
          if (is_curve_bounded(cv,b)) return b;
          X_target_unbounded_curve r;
          cv.current(r);
          const Point& p=r.source(); 
				// Ray is always X_target_unbounded_curve
          out=Bounding_box(
			   Bounding_box(b[0],p+Vector(-num,-num,den))[0],
			   Bounding_box(b[2],p+Vector(num,num,den))[2]);
	  normalize_coordinates(out);
          break;
        }
      case X_curve::LINE:
        {
          if (is_curve_bounded(cv,b)) return b;
          X_unbounded_curve l;
          cv.current(l);
          const Point& p=l.point();
          out=Bounding_box(
            Bounding_box(b[0],p+Vector(-num,-num,den))[0],
            Bounding_box(b[2],p+Vector(num,num,den))[2]);
	  normalize_coordinates(out);
	  break;
        }
      }	
    CGAL_assertion_msg(	
       cv.current_state()==X_curve::SEGMENT||
       cv.current_state()==X_curve::RAY||
       cv.current_state()==X_curve::LINE,
       "Wrong curve type in\nconst Bounding_box\
                           get_bounding_box( \
			   const X_curve& cv,const Bounding_box& b) const");
    CGAL_postcondition(!is_totally_unbounded(out));
#ifdef CGAL_PMBB_DEBUG
      debug_invariant(is_curve_bounded(cv,out));
      debug_invariant(is_bounding_box_bounded(b,out));
#endif
    return out;
  }

  const Bounding_box get_bounding_box(const X_curve& cv,
				      const Ray& ray) const
  /* Returns a bounding box that encloses the intersection of the curve 
     with the ray */
  {
    Point p;
    const Bounding_box& bbox = ref();
    switch(cv.current_state())
      {
      case X_curve::EMPTY:
      case X_curve::POINT: // currently not supported
        return bbox;
      case X_curve::SEGMENT:
        {
          X_bounded_curve s;
          cv.current(s);
          Ray_X_bounded_curve rs(&ray,&s);
          if (!rs.intersection(p)) return bbox;
          break;
        }
      case X_curve::RAY:
        {
          X_source_unbounded_curve r;
          cv.current(r);
          Ray_X_source_unbounded_curve rr(&ray,&r);
          if (!rr.intersection(p)) return bbox;
          break;
        }
      case X_curve::LINE:
        {
          X_unbounded_curve l;
          cv.current(l);
          Ray_X_unbounded_curve rl(&ray,&l);
          if (!rl.intersection(p)) return bbox;
          break;
        }
      }	
    return get_bounding_box(p);
  }

  const Bounding_box get_bounding_box(const X_curve& cv,
				      const Ray& ray,
				      const Bounding_box& b) const
  /* Returns a bounding box that encloses both the curve and the box */
  {
    CGAL_precondition_msg(!is_totally_unbounded(b),"Bounding_box is undefined");
    Point p;
    const Bounding_box& bbox = ref();
    switch(cv.current_state())
      {
      case X_curve::EMPTY:
      case X_curve::POINT: // currently not supported
        return bbox;
      case X_curve::SEGMENT:
        {
          X_bounded_curve s;
          cv.current(s);
          Ray_X_bounded_curve rs(&ray,&s);
          if (!rs.intersection(p)) return bbox;
          break;
        }
      case X_curve::RAY:
        {
          X_source_unbounded_curve r;
          cv.current(r);
          Ray_X_source_unbounded_curve rr(&ray,&r);
          if (!rr.intersection(p)) return bbox;
          break;
        }
      case X_curve::LINE:
        {
          X_unbounded_curve l;
          cv.current(l);
          Ray_X_unbounded_curve rl(&ray,&l);
          if (!rl.intersection(p)) return bbox;
          break;
        }
      }	
    Bounding_box out=Bounding_box(
                       Bounding_box(b[0],p+Vector(-num,-num,den))[0],
                       Bounding_box(b[2],p+Vector(num,num,den))[2]);
    normalize_coordinates(out);
    CGAL_postcondition_msg(!is_totally_unbounded(out),
                           "Bounding_box is undefined");
    return out;
  }

////////////////////////////////
// modifying member functions
////////////////////////////////

 inline void set_bounding_box(const Bounding_box& b) {

#ifdef CGAL_PMBB_DEBUG
   std::cout << "\nset_bounding_box("; debug();
#endif

    ref()=b;

#ifdef CGAL_PMBB_DEBUG
    std::cout << ")->("; debug(); std::cout << ")" << std::flush;
#endif

  }

  inline void set_bounding_box() {
    ref()=unbounded_box();

#ifdef CGAL_PMBB_DEBUG
    std::cout << "\nset_bounding_box("; debug(); 
    std::cout << ")" << std::flush;
#endif

  }

  const Bounding_box& increase_bounding_box(const Point& p,
                                            const Bounding_box& b)
    /* Sets the inner bounding box to encloses both the box and the point */
  {
    
#ifdef CGAL_PM_DEBUG
    std::cout << "\nincrease_bounding_box(" << p << "|"; debug();
    std::cout << ")->(";
#endif
    
    Bounding_box& bbox = ref(); 
    bbox=get_bounding_box(p,b);
    
#ifdef CGAL_PM_DEBUG
    debug(); std::cout << ")";
#endif
    
    return bbox;
  }

  const Bounding_box& increase_bounding_box(const X_curve& cv,
                                            const Bounding_box& b)
    /* Sets the inner bounding box to a bounding box that encloses the 
       curve */
  {
    Bounding_box& bbox = ref();
    
#ifdef CGAL_PM_DEBUG
    std::cout << "\nincrease_bounding_box(" << cv << "|"; debug();
    std::cout << ")->(";
#endif

    bbox=get_bounding_box(cv,b);

#ifdef CGAL_PM_DEBUG
    debug(); std::cout << ")" << std::endl;
#endif

    return bbox;
    }

  const Bounding_box& increase_bounding_box(const X_curve& cv, 
						 const Ray& ray,
						 const Bounding_box& b)
  /* Sets the inner bounding box to a bounding box that encloses the 
     curve */

    {
#ifdef CGAL_PM_DEBUG
      std::cout << "\nincrease_bounding_box(" << cv << "|" << ray << "|";
    debug();
    std::cout << ")->(";
#endif
    Bounding_box& bbox = ref();
      bbox=get_bounding_box(cv,ray,b);
#ifdef CGAL_PM_DEBUG
      debug();
    std::cout << ")";
#endif
    return bbox;
    }

  const Bounding_box& increase_bounding_box(const Point& p)
  /* Enlarges the inner bounding box to enclose the point. */
    {
#ifdef CGAL_PM_DEBUG
      std::cout << "\nincrease_bounding_box(" << p << "|";
    debug();
    std::cout << ")->(";
#endif
    Bounding_box& bbox = ref();
      if (is_totally_unbounded()) // bbox is undefined
	bbox=get_bounding_box(p); 
      else
        increase_bounding_box(p,bbox);
#ifdef CGAL_PM_DEBUG
      debug();
    std::cout << ")";
#endif
    return bbox;
    }

  const Bounding_box& increase_bounding_box(const X_curve& cv)
  /* Sets the inner bounding box to a bounding box that encloses the curve */
    {
#ifdef CGAL_PM_DEBUG
      std::cout << "\nincrease_bounding_box(" << cv << "|";
      debug();
#endif      
      Bounding_box& bbox = ref();
      if (is_totally_unbounded()) // bbox is undefined
	bbox=get_bounding_box(cv); 
      else 
        increase_bounding_box(cv,ref());
#ifdef CGAL_PM_DEBUG      
      std::cout << ")->(";
      debug();
      std::cout << ")" << std::endl;
#endif
      return bbox;
    }

////////////////////////////////
// data member functions
////////////////////////////////

public:
  bool is_bounding_box_bounded(const Bounding_box& s,const Bounding_box& l) 
    const
/* Returns true in the first bounding box is inside the second bounding box. */
    {
      bool sbounded = !is_totally_unbounded(s),lbounded = !is_totally_unbounded(l);
      if (sbounded && lbounded)
        {
          return is_point_really_bounded(s[0],l) && 
          is_point_really_bounded(s[2],l);
        }
      else if (!sbounded) return !lbounded;
      else return true;
    }

  inline bool is_point_bounded(const Point& p, const Bounding_box& b) const
/* Returns true if the input point is bounded inside the interior 
   of the input Bounding_box. */
  {
    CGAL_precondition_msg(!is_totally_unbounded(b),"Bounding_box is undefined");
    return is_bounding_box_bounded(get_bounding_box(p),b);
  }

  /* Returns whether the point is in the bounding box */
  inline bool is_point_bounded(const Point& p) const
  {
    return is_point_bounded(p,bbox);
  }

  /*
  inline bool is_point_bounded_inside(const Point& p) const
  {
    return is_point_bounded_inside(p,bbox);
  }
  */

  bool is_curve_bounded(const X_curve& cv, const Bounding_box& b) const
  {
    // The bounding box is supposed to be canonical, i.e. its sides are 
    // parralel to the axes.
    
    CGAL_precondition_msg(!is_totally_unbounded(b),"Bounding_box is undefined");
    switch(cv.current_state())
      {
      case X_curve::SEGMENT:
        {
          X_bounded_curve s;
          cv.current(s);
          return is_bounding_box_bounded(get_bounding_box(s),b);
        }
      case X_curve::RAY:
        {
          X_target_unbounded_curve r;
          cv.current(r);
          return is_bounding_box_bounded(get_bounding_box(r.source()),b);
        }
      case X_curve::LINE:
        {
          Bounding_box bin=get_real_bounding_box(b);
          X_bounded_curve o;
          X_unbounded_curve l;
          cv.current(l);
	  Bbox_X_unbounded_curve bl(&bin,&l);
          return (bl.intersection(o));
        }
      default:
#ifdef CGAL_PMBB_DEBUG
        debug_invariant(real_state_invariant(cv));
#endif
        break;
      }
    CGAL_assertion_msg(false,"cv.current_state() out or range");
    return false;
  }

  bool is_curve_bounded(const X_curve& cv) const
  {
    return is_curve_bounded(cv,ref());
  }

  ///////////////////////////////////////////////////

  Segment create_segment(const Point& p,const Point& q) const
  {
    return Segment(p,q);
  }
  void get_corner(const Point& p1,const Point& p2,X_curve& cv1,X_curve& cv2, const Boundary_type& b) const
  {
    Bounding_box box = Bounding_box(p1,p2);
    switch (b)
      {
      case TOP:
        cv1=get_x_curve_boundary(LEFT,box);
        cv2=get_x_curve_boundary(TOP,box);
        break;
      case RIGHT:
        cv1=get_x_curve_boundary(TOP,box);
        cv2=get_x_curve_boundary(RIGHT,box);
        break;
      case BOTTOM:
        cv1=get_x_curve_boundary(RIGHT,box);
        cv2=get_x_curve_boundary(BOTTOM,box);
        break;
      case LEFT:
        cv1=get_x_curve_boundary(BOTTOM,box);
        cv2=get_x_curve_boundary(LEFT,box);
        break;
      }
  }
  /* Postcondition: cv is the disjoint union of cv1 and cv2 */
  void split_curve(const X_curve& cv, X_curve& cv1, X_curve& cv2, 
		   const Point& p) const
    {
      switch(cv.current_state())
	{
	case X_curve::POINT:
	  CGAL_warning(cv.current_state()!=X_curve::POINT);
	  return;
	case X_curve::EMPTY:
	  CGAL_warning(cv.current_state()!=X_curve::EMPTY);
	  return;
	case X_curve::SEGMENT:
	  {
	    X_bounded_curve s;
	    cv.current(s);
#ifdef CGAL_PM_DEBUG
	    CGAL_warning(s.has_on(p));
	    CGAL_warning(s.collinear_has_on(p));
#endif
	    cv1=X_bounded_curve(s.source(),p);
	    cv2=X_bounded_curve(p,s.target());
	    break;
	  }
	case X_curve::RAY:
	  {
	    X_target_unbounded_curve r;
	    cv.current(r);
#ifdef CGAL_PM_DEBUG
	    CGAL_warning(r.has_on(p));
	    CGAL_warning(r.collinear_has_on(p));
#endif
	    cv1=X_bounded_curve(r.source(),p);
	    cv2=X_target_unbounded_curve(p,r.direction());
	    break;
	  }
	case X_curve::LINE:
	  {
	    X_unbounded_curve l;
	    cv.current(l);
#ifdef CGAL_PM_DEBUG
	    CGAL_warning(l.has_on_boundary(p));
#endif
	    const X_target_unbounded_curve 
	      ray(p,l.direction());
	    cv1=X_curve(ray.opposite(),false);
	    cv2=ray;
	    break;
	  }
	}
    }
  
 public:  
  static const Bounding_box unbounded_box();
 private:
  static const Bounding_box unbounded_box_;

#if defined(CGAL_PM_DEBUG) || defined(CGAL_PMBB_DEBUG)
public:
  void debug() const
  {
		debug(ref());    
  }
  void debug(const Bounding_box& b) const
  {
    std::cout << "Bounding_box[";
    if (is_totally_unbounded(b))
      std::cout << "Plane";
    else
      std::cout << b;
    std::cout << "]" << std::flush;
  }
#endif

#ifdef CGAL_PMBB_DEBUG
  bool debug_invariant(bool flag) const{
    CGAL_assertion(flag);
    return flag;
  }
  bool bound_state_invariant(const X_curve& cv) const
  {
    return cv.bound_state()==X_curve::NO_UNBOUNDED||
      cv.bound_state()==X_curve::MIN_UNBOUNDED||
      cv.bound_state()==X_curve::MAX_UNBOUNDED||
      cv.bound_state()==X_curve::BOTH_UNBOUNDED;
  }
  bool real_state_invariant(const X_curve& cv) const
  {
    return cv.bound_state()==X_curve::SEGMENT||
      cv.bound_state()==X_curve::RAY||
      cv.bound_state()==X_curve::LINE;
  }
  bool is_bounded_invariant(const Point& p, const Bounding_box& b) const
  {
    return is_point_bounded(p,b);
  }
  bool is_bounded_invariant(const X_curve& cv, const Bounding_box& b) const
  {
    return is_curve_bounded(cv,b);
  }
#endif

private:
  inline const Bounding_box& ref() const {return ref(bbox);}
  inline const Bounding_box& ref(const Indirect_bounding_box& bbox) const 
  {return bbox.pointer()->e0;}
  inline Bounding_box& ref() {return ref(bbox);}
  inline Bounding_box& ref(const Indirect_bounding_box& bbox)
  {return bbox.pointer()->e0;}
  inline bool is_totally_unbounded() const {return is_totally_unbounded(ref());}
  inline bool is_totally_unbounded(const Bounding_box& b) const
  {
    return  b.identical(unbounded_box());
  }

protected:
  /*
  bool is_point_really_on_boundary(const Point& p, const Bounding_box& b) const
    {
      CGAL_precondition_msg(!is_totally_unbounded(b),
			    "Bounding_box is undefined");
      return b.has_on_boundary(p);
    }
  inline bool is_point_really_on_boundary(const Point& p) const
// Returns whether the point is on the bounding box inside
  {
    return is_point_on_boundary(p,ref());
  }
  */
  inline bool is_point_really_bounded(const Point& p, const Bounding_box& b) const
    /* Returns true if the input point is inside the closed input 
       Bounding_box. */
  {
    CGAL_precondition_msg(!is_totally_unbounded(b),"Bounding_box is undefined");
    return b.bounded_side(p)!=ON_UNBOUNDED_SIDE;
  }
  Bounding_box get_real_bounding_box(const Bounding_box& b) const
    /* used internally to calculate the interior of a bounding box. */
  {
    return Bounding_box(b[0]+Vector(num,num,den),b[2]+Vector(-num,-num,den));
  }

protected:
  Indirect_bounding_box bbox;
  RT num,den;
  
};

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION 
#include <CGAL/Pm_straight_exact_traits.C>
#endif

#endif // CGAL_PM_STRAIGHT_EXACT_TRAITS_H



