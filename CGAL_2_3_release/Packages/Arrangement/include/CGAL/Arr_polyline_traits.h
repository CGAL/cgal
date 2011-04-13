// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/Arr_polyline_traits.h
// package       : Arrangement (1.77)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Iddo Hanniel
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_POLYLINE_TRAITS_H
#define CGAL_ARR_POLYLINE_TRAITS_H

#include <CGAL/basic.h>

#include <list>
#include <deque>
#include <vector>

#include <algorithm>

#include <CGAL/predicates_on_points_2.h>
#include <CGAL/predicates_on_lines_2.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>

#include <CGAL/Point_2.h>
#include <CGAL/squared_distance_2.h>

#include <CGAL/Pm_segment_exact_traits.h>

//#include <typeinfo>

CGAL_BEGIN_NAMESPACE


template <class R,
  class Container  = std::vector<Point_2<R> >
>
class Arr_polyline_traits {
public:
  typedef Arr_polyline_traits<R> Self;

  typedef Point_2<R> Point;
  typedef Vector_2<R> Vector; //for drivative

  typedef Container Curve;
  typedef Container X_curve;

  typedef enum
  {
    UNDER_CURVE = -1,
    ABOVE_CURVE = 1,
    ON_CURVE = 2,
    CURVE_NOT_IN_RANGE = 0
    //CURVE_VERTICAL = 3
  } Curve_point_status;


  Arr_polyline_traits() {
  }

  Comparison_result compare_x(const Point& p0, const Point& p1) const {
    return CGAL::compare_x(p0,p1);
  }
  Comparison_result compare_y(const Point& p0, const Point& p1) const {
    return CGAL::compare_y(p0,p1);
  }

  //on X_curve only - not Curve!
  bool curve_is_vertical(const X_curve& cv) const {
    CGAL_assertion(is_x_monotone(cv));
    return compare_x(curve_source(cv),curve_target(cv))==EQUAL;
  } 

  bool curve_is_in_x_range(const X_curve& cv, const Point& p) const {
    CGAL_assertion(is_x_monotone(cv));
    return (compare_x(p,curve_source(cv)) * 
	    compare_x(p,curve_target(cv))) <= 0 ;
  }

  Comparison_result curve_compare_at_x(const X_curve& cv1, 
				       const X_curve& cv2, 
				       const Point& p) const {
    CGAL_assertion(is_x_monotone(cv1));
    CGAL_assertion(is_x_monotone(cv2));

    if (!curve_is_in_x_range(cv1,p) || !curve_is_in_x_range(cv2,p) )
      return EQUAL;

    if (curve_is_vertical(cv1) && curve_is_vertical(cv2))
      return EQUAL; //otherwise - compare_y_at_x throws an assertion 

    typename X_curve::const_iterator pit_1   = cv1.begin(),
      pit_2   = cv2.begin();
    typename X_curve::const_iterator after_1 = pit_1,
      after_2 = pit_2;
    ++after_1; ++after_2;

    // read: as long as both *pit and *after are on the same
    //       side of p 
    for( ; (compare_x(*pit_1,p) * compare_x(*after_1,p)) /*<=*/ >0 ; 
	 ++pit_1,++after_1) {}
    for( ; (compare_x(*pit_2,p) * compare_x(*after_2,p)) > 0 ; 
	 ++pit_2,++after_2) {}

    // the R here is the template parameter (see class def. above)
    Pm_segment_exact_traits<R> segment_traits;
    
    const typename Pm_segment_exact_traits<R>::X_curve
      seg1(*pit_1, *after_1),
      seg2(*pit_2, *after_2);

    // bug fix, shai 19/03/2000:
    // the polylines may contain vertical segments
    return (segment_traits.curve_compare_at_x(seg1, seg2, p));
  }

  //precondition - x-monotone
  Comparison_result curve_compare_at_x_left(const X_curve& cv1, 
					    const X_curve& cv2,
					    const Point& p) const {
    CGAL_assertion(is_x_monotone(cv1));
    CGAL_assertion(is_x_monotone(cv2));


    if (!curve_is_in_x_range(cv1,p) || !curve_is_in_x_range(cv2,p) )
      return EQUAL;

    Point leftmost1=(compare_x(curve_source(cv1),curve_target(cv1))==LARGER) ?
      curve_target(cv1) : curve_source(cv1);
    Point leftmost2=(compare_x(curve_source(cv2),curve_target(cv2))==LARGER) ?
      curve_target(cv2) : curve_source(cv2);

    //special cases wher returns EQUAL
    if (curve_is_vertical(cv1) || (curve_is_vertical(cv2))) return EQUAL;
    if (compare_x(leftmost1,p)!=SMALLER || compare_x(leftmost2,p)!=SMALLER) {
      return EQUAL;
    }



    
    typename X_curve::const_iterator pit=cv1.begin();
    typename X_curve::const_iterator after=pit; ++after;
    
    for( ; (compare_x(*pit,p) * compare_x(*after,p)) > 0 ; ++pit,++after) {}
    
    Line_2<R> l1(*pit,*after);
    
    pit=cv2.begin();
    after=pit; ++after;
    for( ; (compare_x(*pit,p) * compare_x(*after,p)) > 0 ; ++pit,++after) {}
    
    Line_2<R> l2(*pit,*after);
    
    Comparison_result r=CGAL::compare_y_at_x(p,l1,l2); 
    
    if ( r != EQUAL)
      return r;
    else {
      // check if they are right endpoints (and compare to the one from 
      // the left) otherwise -
      return (CGAL::compare_y_at_x(point_to_left(p),l1,l2)); 
    }

  } 



  Comparison_result curve_compare_at_x_right(const X_curve& cv1, 
					     const X_curve& cv2,
					     const Point& p) const {
    CGAL_assertion(is_x_monotone(cv1));
    CGAL_assertion(is_x_monotone(cv2));

    if (!curve_is_in_x_range(cv1,p) || !curve_is_in_x_range(cv2,p))
      return EQUAL;

    Point rightmost1=(compare_x(curve_source(cv1),curve_target(cv1))==SMALLER)
      ? curve_target(cv1) : curve_source(cv1);
    Point rightmost2=(compare_x(curve_source(cv2),curve_target(cv2))==SMALLER)
      ? curve_target(cv2) : curve_source(cv2);

    //special cases wher returns EQUAL
    if (curve_is_vertical(cv1) || (curve_is_vertical(cv2))) return EQUAL;
    if (compare_x(rightmost1,p)!=LARGER || compare_x(rightmost2,p)!=LARGER) {
      return EQUAL;
    }

    //not defined right of curve:
    //    CGAL_assertion(!is_same(rightmost1,p) && !is_same(rightmost2,p)); 


    if (!curve_is_in_x_range(cv1,p) || !curve_is_in_x_range(cv2,p) )
      return EQUAL;

    typename X_curve::const_iterator pit=cv1.begin();
    typename X_curve::const_iterator after=pit; ++after;
    
    for( ; (compare_x(*pit,p) * compare_x(*after,p)) > 0 ; ++pit,++after) {}
    
    Line_2<R> l1(*pit,*after);
    
    pit=cv2.begin();
    after=pit; ++after;
    for( ; (compare_x(*pit,p) * compare_x(*after,p)) > 0 ; ++pit,++after) {}

    Line_2<R> l2(*pit,*after);
    
    Comparison_result r=CGAL::compare_y_at_x(p,l1,l2); 
    
    if ( r != EQUAL) {
      return r;
    }
    else
      // check if they are left endpoints (and compare to the one from the 
      // right) otherwise -
      
      //debug
      { 
	CGAL_assertion(curve_compare_at_x(cv1,cv2,p) == EQUAL);
	return (CGAL::compare_y_at_x(point_to_right(p),l1,l2)); 
      }
    //debug
    
  }
  


  Curve_point_status 
  curve_get_point_status(const X_curve &cv, const Point& p) const
  {
    CGAL_assertion(is_x_monotone(cv));

    if (!curve_is_in_x_range(cv, p))
      return CURVE_NOT_IN_RANGE;
    if (curve_is_vertical(cv)) {
      if (compare_y(curve_source(cv),p)*compare_y(curve_target(cv),p)<=0)
        return ON_CURVE;
      if (compare_y(curve_source(cv),p)==LARGER)
	//bug fix (2/11)
	//return ABOVE_CURVE;
	return UNDER_CURVE;
      if (compare_y(curve_source(cv),p)==SMALLER)
	//return UNDER_CURVE;
	return ABOVE_CURVE;
    }

	typename X_curve::const_iterator pit=cv.begin(),after=pit; ++after;
    while ( (compare_x(*pit,p) * compare_x(*after,p)) > 0 ) {
      ++pit; ++after;
    }
 
    Line_2<R> l(*pit,*after);

    Comparison_result res = CGAL::compare_y_at_x(p, l);

    if (res == SMALLER)
      return UNDER_CURVE;
    if (res == LARGER)
      return ABOVE_CURVE;
    if (res == EQUAL)
      return ON_CURVE;
    return ON_CURVE;
  }
  

  //precondition - same as in pm
  bool curve_is_between_cw(const X_curve& cv,const X_curve& first,
                           const X_curve& second, const Point& p) const
  {
    CGAL_assertion(is_x_monotone(cv));
    CGAL_assertion(is_x_monotone(first));
    CGAL_assertion(is_x_monotone(second));

    CGAL_precondition(is_same(curve_source(cv),p) || 
		      is_same(curve_target(cv),p));
    CGAL_precondition(is_same(curve_source(first),p) || 
		      is_same(curve_target(first),p));
    CGAL_precondition(is_same(curve_source(second),p) || 
		      is_same(curve_target(second),p));



    X_curve cv0 = first;
    X_curve cv1 = second;
    X_curve cvx = cv;



    if ( !is_same(curve_source(cv0),p) ) cv0 = curve_flip(cv0);
    if ( !is_same(curve_source(cv1),p) ) cv1 = curve_flip(cv1);
    if ( !is_same(curve_source(cvx),p) ) cvx = curve_flip(cvx);

    
    typename X_curve::iterator xcit=cv0.begin();++xcit;
    Point p0(*xcit);
    xcit=cv1.begin(); ++xcit;
    Point p1(*xcit);
    xcit=cvx.begin(); ++xcit;
    Point px(*xcit);

    if (is_same(p0,p1))
      return true; 
    
    int or0=orientation(p0,p,px);
    int or1=orientation(p1,p,px);
    // Bug Fix, Shai, Jan, 8, 2001
    // 'or' is a keyword in C++, changed to 'orient'
    int orient=or0*or1;
    
    if (orient < 0) { //one is a leftturn the other rightturn
      return (or0 == LEFTTURN); //leftturn
    }
    else { //both are either left or right turns (or one is colinear)
      return (orientation(p0,p,p1)==RIGHTTURN); //rightturn
    }
  }
  


  bool curve_is_same(const X_curve& cv1, const X_curve& cv2) const {
    CGAL_assertion(is_x_monotone(cv1));
    CGAL_assertion(is_x_monotone(cv2));

    if (cv1.size()!=cv2.size())
      return false;
    typename X_curve::const_iterator it1=cv1.begin();
    typename X_curve::const_iterator it2=cv2.begin();
    while (it1!=cv1.end()) {
      if (!is_same(*it1++,*it2++))
        return false;
    }
    return true;
  }

  // SHAI: shouldn't cv  be of type Curve (even if X_curve is the same)
  Point curve_source(const X_curve& cv) const {
    return *(cv.begin());
  }

  // SHAI: shouldn't cv  be of type Curve (even if X_curve is the same)
  Point curve_target(const X_curve& cv) const {
    //debug
    //this seems to not work with vector
    //return *(--cv.end());
    typename X_curve::const_iterator it=cv.end(); --it;
    return *it;
  }

  Point point_to_left(const Point& p) const {
    return p+Vector(-1,0);;
  }
  Point point_to_right(const Point& p) const {
    return p+Vector(1,0);;
  }
  

  ///////////////////////////////////////////////////////
  //         ARRANGEMENT FUNCS



  X_curve curve_flip(const X_curve& cv) const {
    X_curve cv1(cv);
    //cv1.reverse();
    std::reverse(cv1.begin(),cv1.end()); 
    return cv1; 
  }

  bool is_x_monotone(const Curve& cv) const {
    CGAL_assertion(cv.size()>=2); //one point is not a curve

    if (cv.size()==2) return true; //segment

    typename X_curve::const_iterator p0=cv.begin();
    typename X_curve::const_iterator p1=p0; ++p1;
    typename X_curve::const_iterator p2=p1; ++p2;

    for(; p2!=cv.end(); ++p0,++p1,++p2) {
      if ( compare_x(*p0,*p1) * compare_x(*p1,*p2) <=0 )
        return false;
    }

    // <= a matter of decision - only one vertical segment is considered 
    // x-monotone
    return true;

  }

  //cuts into x-monotone curves, each vertical segment is 1 x-monotone curve
  //and not part of a bigger x-monotone polyline
  void make_x_monotone(const Curve& cv, std::list<Curve>& l) {
    CGAL_assertion(cv.size()>=2); //one point is not a curve

    if (cv.size()==2) { //segment
      l.push_back(cv);
      return;
    }

    typename X_curve::const_iterator p0=cv.begin();
    typename X_curve::const_iterator p1=p0; ++p1;
    typename X_curve::const_iterator p2=p1; ++p2;

    typename X_curve::const_iterator last_cut=p0;

    for(; p2!=cv.end(); ++p0,++p1,++p2) {
      //in future use constants instead of compare_x
      if (compare_x(*p0,*p1)==EQUAL) { //vertical segment - (p0,p1) 
        if (p0!=last_cut) { //needed for the case:
          //    | p1
          //    o p0=lastcut (was p1 before)
          //    |

          l.push_back(X_curve(last_cut,p1)); //constructor deque(first,beyond)
          //push_back the curve (last_cut...p0)
        }
        l.push_back(X_curve(p0,p2)); //push_back the segment (p0,p1)
        last_cut=p1;
      }
      else
        if ( compare_x(*p0,*p1) * compare_x(*p1,*p2) <= 0 ) {
          l.push_back(X_curve(last_cut,p2));
          last_cut=p1;
        }
    }

    l.push_back(X_curve(last_cut,p2)); //push the residue (last cut to end)


    CGAL_assertion(p2==cv.end());
    
  }


  void curve_split(const X_curve& cv, X_curve& c1, X_curve& c2, 
		   const Point& split_pt) {
    
    //split curve at split point into c1 and c2
    CGAL_precondition(curve_get_point_status(cv,split_pt)==ON_CURVE);
    CGAL_precondition(CGAL::compare_lexicographically_xy(curve_source(cv),
							 split_pt) != EQUAL);
    CGAL_precondition(CGAL::compare_lexicographically_xy(curve_target(cv),
							 split_pt) != EQUAL);

    typename X_curve::const_iterator p0=cv.begin();
    typename X_curve::const_iterator p1=p0; ++p1;
 
    bool split_at_vertex=false;

    for (; p1 != cv.end(); ++p0,++p1) {
      if (is_same(*p1,split_pt)) {
	split_at_vertex=true;
	break;
      }
 
      if (compare_x(*p0,split_pt) * compare_x(split_pt,*p1) >= 0) {
	//in x range 
        break;
      }
    } 
    
    c1.clear(); c2.clear();

    typename X_curve::const_iterator ci=cv.begin();
    while (ci!=p1) {
      c1.push_back(*ci);
      ++ci;
    }

    if (!split_at_vertex) {
      c1.push_back(split_pt);
      c2.push_back(split_pt);
    }
    else {
      c1.push_back(*p1);
    }

    while (ci!=cv.end()) {
      c2.push_back(*ci);
      ++ci;
    }

    //moved this up to enable use of vector as container
    /*
      if (!split_at_vertex) {
      c1.push_back(split_pt);
      c2.push_front(split_pt);
      }
      else {
      c1.push_back(*p1);
      }
    */


  }

public:

  //returns true iff the intersectionis lexicographically strictly right of pt

  bool do_intersect_to_right(const X_curve& ca, const X_curve& cb,
                             const Point& pt) const 
  {
    CGAL_assertion(is_x_monotone(ca));
    CGAL_assertion(is_x_monotone(cb));


    // check if both first points are left of pt, if they are reach the 
    // points directly left of pt, and check if their segments intersect  
    //to the right of pt, if not continue with a normal sweep until 
    //we find an intersection point or we reach the end.
    
    //do a flip or can we assume they are left to right ??
    X_curve c1(ca),c2(cb);
    if (lexicographically_xy_larger(curve_source(ca), curve_target(ca)) ==
	LARGER )
      c1=curve_flip(ca);
    if (lexicographically_xy_larger(curve_source(cb),curve_target(cb)) == 
	LARGER )
      c2=curve_flip(cb);

    typename X_curve::const_iterator i1s=c1.begin(),i1e=c1.end();
    typename X_curve::const_iterator i1t=i1s; ++i1t;

    typename X_curve::const_iterator i2s=c2.begin(),i2e=c2.end();
    typename X_curve::const_iterator i2t=i2s; ++i2t;

    int number_to_left=0; //increment this variable if curve starts left of pt

    if (!lexicographically_xy_larger (*i1s,pt)) {

      //increment to nearest from the left of pt
      ++number_to_left;
      for (; i1t!=i1e; ++i1s,++i1t) {
	if (lexicographically_xy_larger (*i1t,pt)) break;
      }
      if (i1t==i1e) return false; //c1 ends to the left of pt
    }
    //now i1s holds the source vertex and i1t holds the target

    if (!lexicographically_xy_larger (*i2s,pt)) {
      //increment 
      ++number_to_left;
      for (; i2t!=i2e; ++i2s,++i2t) {
	if (lexicographically_xy_larger (*i2t,pt)) break;
      }
      if (i2t==i2e) return false; //c2 ends to the left of pt
    }

    if (number_to_left==2) {
      //check if intersection exists and is lex larger
      Object result;
      Point i_pt;
      Segment_2<R> i_seg;
      
      result = intersection(Segment_2<R>(*i1s,*i1t),
			    Segment_2<R>(*i2s,*i2t));
      if (assign(i_pt,result)) {
        //check if intersection point to the right of pt
        if (lexicographically_xy_larger (i_pt,pt)) 
	  {
	    return true;
	  }
      }
      else

	if (assign(i_seg,result)) {
	  //check if intersection seg to the right of pt
	  if (lexicographically_xy_larger (i_seg.source(),pt) ||
	      lexicographically_xy_larger (i_seg.target(),pt))
	    {
	      return true;
	    }
	}
      //debug
	else {
	  //cerr << "segments don't intersect ??" << endl;
	}
      //advance to the nearer point
      if (lexicographically_xy_larger (*i2t,*i1t)) {
	++i1s; ++i1t;
        if (i1t==i1e) return false;
      }
      else {
	++i2s; ++i2t;
        if (i2t==i2e) return false;
      }

    }
    //NOW we can start sweeping the chains

    while (1) {
      //check for intersection of the segments
      if (do_intersect(Segment_2<R>(*i1s,*i1t),
		       Segment_2<R>(*i2s,*i2t)))
	{
	  return true;
	}
      
      //advance to the nearer point
      if (lexicographically_xy_larger (*i2t,*i1t)) {
	++i1s; ++i1t;
        if (i1t==i1e) return false;
      }
      else {
	++i2s; ++i2t;
        if (i2t==i2e) return false;
      }
    }

  }


  //NOTE: when there is an overlap we will always return a SEGMENT (i.e.,
  //      p1 and p2 will be on a segment) even if the overlap is a polyline
  //      , this is still sufficient for the arrangement. might be
  //      changed in the future.
  bool nearest_intersection_to_right(const X_curve& cv1,
                                     const X_curve& cv2,
                                     const Point& pt,
                                     Point& p1,
                                     Point& p2) const 
  {      
    CGAL_assertion(is_x_monotone(cv1));
    CGAL_assertion(is_x_monotone(cv2));

    bool found( false);

    // bug fix:
    // curves do not necessarily intersect
    if ( ! do_intersect_to_right(cv1,cv2,pt)) return false;

    X_curve c1(cv1),c2(cv2);
    if ( ! lexicographically_xy_smaller (curve_source(c1),curve_target(c1)))
      c1=curve_flip(cv1);
    if ( ! lexicographically_xy_smaller (curve_source(c2),curve_target(c2)))
      c2=curve_flip(cv2);
      
    // check if both first points are left of pt, if they are reach the 
    // points directly left of pt, and check if their segments intersect  
    //to the right of pt, if not continue with a normal sweep until 
    //we find an intersection point or we reach the end.
      
    typename X_curve::const_iterator i1s=c1.begin(),i1e=c1.end();
    typename X_curve::const_iterator i1t=i1s; ++i1t;
      
    typename X_curve::const_iterator i2s=c2.begin(),i2e=c2.end();
    typename X_curve::const_iterator i2t=i2s; ++i2t;
      
    int number_to_left=0; //increment this variable if curve starts left of pt
      
    if (!lexicographically_xy_larger (*i1s,pt)) {
      //increment to nearest from the left of pt
      ++number_to_left;
      for (; i1t!=i1e; ++i1s,++i1t) {
	if (lexicographically_xy_larger (*i1t,pt)) break;
      }
      if (i1t==i1e) return false;
    }
      
    //now i1s holds the source vertex and i1t holds the target
      
    if (!lexicographically_xy_larger (*i2s,pt)) {
      //increment 
      ++number_to_left;
      for (; i2t!=i2e; ++i2s,++i2t) {
	if (lexicographically_xy_larger (*i2t,pt)) break;
      }
      if (i2t==i2e) return false ; //c2 ends to the left of pt
    }
      
    if (number_to_left==2) {
      //check if intersection exists and is lex larger
      Object result=intersection(Segment_2<R>(*i1s,*i1t),
				 Segment_2<R>(*i2s,*i2t));
      Point i_pt;
      Segment_2<R> i_seg;
      if (assign(i_pt,result)) {
	//check if intersection point to the right of pt
	if (lexicographically_xy_larger (i_pt,pt)) {
	  //debug
#ifndef ARR_USES_LEDA_RATIONAL  //normalize if we are with rational numbers
	  p1=p2=i_pt;
	  found =  true;
#else
	  p1=p2=Point(i_pt.x().normalize(),i_pt.y().normalize());
	  found = true;
#endif
	}
      }
       
      if ( ! found && assign(i_seg,result)) {
	//check if intersection seg to the right of pt
	if (lexicographically_xy_larger (i_seg.source(),i_seg.target()))
	  i_seg=i_seg.opposite();

	if (lexicographically_xy_larger (i_seg.target(),pt)) {
	  p2=i_seg.target();
	  if (lexicographically_xy_larger (i_seg.source(),pt)) {
	    p1=i_seg.source();}
	  else {
	    // p1=pt;
            // Modified by Eug
            // Performing vertical ray shooting from pt.
	    // Finding the intersection point. We know by now
	    // that there is exactly ONE point. Assinging this
	    // point to p1.
            Point ap1( pt.x(), i_seg.source().y() );
            Point ap2( pt.x(), i_seg.target().y() );
            Segment_2<R> vertical_pt_x_base( ap1, ap2 );
            Object i_obj = intersection( vertical_pt_x_base, i_seg );
            assign( p1, i_obj );
          }         
	  found = true;
	}
      }
        
      if ( ! found) {
	//advance to the nearer point
	if (lexicographically_xy_larger (*i2t,*i1t)) {
	  ++i1s; ++i1t;
	  if (i1t==i1e) return false;
	}
	else {
	  ++i2s; ++i2t;
	  if (i2t==i2e) return false;
	}
      }
        
    }
      
    //NOW we can start sweeping the chains
      
    while ( ! found) {
      //check intersection of the segments
        
      Object result;
        
      result = intersection(Segment_2<R>(*i1s,*i1t),
			    Segment_2<R>(*i2s,*i2t));
        
      Point i_pt;
      Segment_2<R> i_seg;
      if (assign(i_pt,result)) {
          
#ifndef ARR_USES_LEDA_RATIONAL  //normalize if we are with rational numbers
	p1=p2=i_pt;
#else
	p1=p2=Point(i_pt.x().normalize(),i_pt.y().normalize());
#endif
	found = true;
      }
        
       
      if (!found && assign(i_seg,result)) {
	//check if intersection seg to the right of pt
	if (lexicographically_xy_larger (i_seg.source(),i_seg.target()))
          i_seg=i_seg.opposite();

	if (lexicographically_xy_larger (i_seg.target(),pt)) {
	  p2=i_seg.target();
	  if (lexicographically_xy_larger (i_seg.source(),pt)) {
	    p1=i_seg.source();}
	  else {
	    // p1=pt;
            // Modified by Eug
            // Performing vertical ray shooting from pt.
	    // Finding the intersection point. We know by now
	    // that there is exactly ONE point. Assinging this
	    // point to p1.
            Point ap1( pt.x(), i_seg.source().y() );
            Point ap2( pt.x(), i_seg.target().y() );
            Segment_2<R> vertical_pt_x_base( ap1, ap2 );
            Object i_obj = intersection( vertical_pt_x_base, i_seg );
            assign( p1, i_obj );
	  }    
	  found = true;
	}
      }
        

      if ( ! found) {
	//advance to the nearer point
	if (lexicographically_xy_larger (*i2t,*i1t)) {
	  ++i1s; ++i1t;
	  if (i1t==i1e) return false;
	}
	else {
	  ++i2s; ++i2t;
	  if (i2t==i2e) return false;
	}
      }
        
    } // end while ( ! found)
    
    // if the x point is at the end of a segment, then there might be 
    // an overlap in the continious of the polyline
    if ( found && is_same( p1, p2) ) {
      typename X_curve::const_iterator s1=i1s, t1=i1t, s2=i2s, t2=i2t;
      s1++; t1++; s2++; t2++;
      if (t1 != i1e && t2 != i2e) {

	// check for overlap after x point
	Object result;
	Segment_2<R> i_seg;

	if ( is_same( p1, *i1t) && is_same( p1, *i2t)) {
	  result = intersection(Segment_2<R>(*s1, *t1), 
				Segment_2<R>(*s2, *t2));
	}
	else if ( is_same( p1, *i1t)) {
	  result = intersection(Segment_2<R>(*s1, *t1),
				Segment_2<R>(*i2s, *i2t));
	}
	else if ( is_same( p1, *i2t)) {
	  result = intersection(Segment_2<R>(*i1s, *i1t),
				Segment_2<R>(*s2, *t2));
	}
	if (assign(i_seg,result)) {
	  // no need to check whether intersection seg to the right of pt, 
	  // because we have already found an x point to the right of pt.
	  if (lexicographically_xy_larger (i_seg.source(),i_seg.target()))
	    i_seg=i_seg.opposite();
	  p1 = i_seg.source();
	  p2 = i_seg.target();
	}
      }  
    }

    return found;
  }


  bool curves_overlap(const X_curve& ca, const X_curve& cb) const {
    CGAL_assertion(is_x_monotone(ca));
    CGAL_assertion(is_x_monotone(cb));

    //do a flip so they are left to right 
    X_curve c1(ca),c2(cb);
    if (lexicographically_xy_larger(curve_source(ca),curve_target(ca)))
      c1=curve_flip(ca);
    if (lexicographically_xy_larger(curve_source(cb),curve_target(cb)))
      c2=curve_flip(cb);

    typename X_curve::const_iterator i1s=c1.begin(),i1e=c1.end();
    typename X_curve::const_iterator i1t=i1s; ++i1t;

    typename X_curve::const_iterator i2s=c2.begin(),i2e=c2.end();
    typename X_curve::const_iterator i2t=i2s; ++i2t;

    //now i1s holds the source vertex and i1t holds the target
    Point i_pt;
    Segment_2<R> i_seg;
    Segment_2<R> s1(*i1s,*i1t),s2(*i2s,*i2t);
    Object res=intersection(s1,s2);
    
    if (assign(i_seg,res)) {
      if (!is_same(i_seg.source(),i_seg.target()))
	return true;
    }
    //advance to the nearer point
    if (lexicographically_xy_larger(*i2t,*i1t)) {
      ++i1s; ++i1t;
      if (i1t==i1e) return false;
    }
    else {
      ++i2s; ++i2t;
      if (i2t==i2e) return false;
    }

    //NOW we can start sweeping the chains

    while (1) {
      Segment_2<R> i_seg;
      Segment_2<R> s1(*i1s,*i1t),s2(*i2s,*i2t);
      res=intersection(s1,s2);
      if (assign(i_seg,res)) {
	if (!is_same(i_seg.source(),i_seg.target()))
	  return true;
      }

      if (lexicographically_xy_larger(*i2t,*i1t)) {
        ++i1s; ++i1t;
        if (i1t==i1e) return false;
      }
      else {
        ++i2s; ++i2t;
        if (i2t==i2e) return false;
      }
    }    
  }


  X_curve curve_reflect_in_x_and_y( const X_curve& cv) const
  {
    X_curve reflected_cv;
    typename Curve::const_iterator it  = cv.begin();
    for  (; it != cv.end(); it++)
      {
        reflected_cv.push_back( point_reflect_in_x_and_y( *it)); 
      }
    return reflected_cv;
  }


  Point point_reflect_in_x_and_y (const Point& pt) const
  {
    // use hx(), hy(), hw() in order to support both Homogeneous and Cartesian
    Point reflected_pt( -pt.hx(), -pt.hy(), pt.hw());
    return reflected_pt;
  }


  ////////////////////////////////////////////////////////////////////
  // PRIVATE
private:


  bool is_same(const Point &p1, const Point &p2) const
  {
    return (compare_x(p1, p2) == EQUAL) &&
      (compare_y(p1, p2) == EQUAL);
  }


public:
  void display(const X_curve& cv) const
  {
    typename X_curve::const_iterator it=cv.begin(),eit=cv.end();
    while(it!=eit) { std::cerr << *it++;}
  }
  
  //the same for window stream
};


CGAL_END_NAMESPACE




#endif





