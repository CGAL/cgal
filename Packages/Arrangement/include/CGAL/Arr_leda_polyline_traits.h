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
// release       : $CGAL_Revision: CGAL-2.4-I-63 $
// release_date  : $CGAL_Date: 2002/03/15 $
//
// file          : include/CGAL/Arr_leda_polyline_traits.h
// package       : Arrangement (2.37)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Iddo Hanniel
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_ARR_LEDA_POLYLINE_TRAITS_H
#define CGAL_ARR_LEDA_POLYLINE_TRAITS_H

#include <CGAL/LEDA_basic.h>
#include <CGAL/basic.h>

#include <list>
#include <deque>
#include <vector>

#include <algorithm>
#include <CGAL/rat_leda_in_CGAL_2.h>
#include <LEDA/rat_point.h>

#include <CGAL/Pm_segment_traits_leda_kernel_2.h>

//the following is for a type check (creates compiler problems,not implemented)
//#include <typeinfo>

CGAL_BEGIN_NAMESPACE


template <class FT_, class Container  = std::vector<leda_rat_point > >
class Arr_leda_polyline_traits : public Pm_segment_traits_leda_kernel_2<FT_>
{
public:
  typedef FT_                                           FT;
  typedef Arr_leda_polyline_traits<FT,Container>        Self;

  typedef leda_rat_point                                Point_2;
  typedef leda_rat_vector                               Vector_2;
  typedef leda_rat_segment                              Segment_2;

  typedef Container                                     Curve_2;
  typedef Container                                     X_curve_2;

  // Obsolete, for backward compatibility
  typedef Point_2                                       Point;
  typedef Vector_2                                      Vector;
  typedef X_curve_2                                     X_curve;
  typedef Curve_2                                       Curve;
  typedef Segment_2                                     Segment;

  typedef enum
  {
    UNDER_CURVE = -1,
    ABOVE_CURVE = 1,
    ON_CURVE = 2,
    CURVE_NOT_IN_RANGE = 0
  } Curve_point_status;


  Arr_leda_polyline_traits() { }

  Comparison_result compare_x(const Point_2& p0, const Point_2& p1) const {
    return (Comparison_result)(

#if (__LEDA__ >= 380)
                               // full specification is required for msvc:
                               Self::Point_2::cmp_x(p0,p1)
#else // backward compatability to LEDA   
                               compare(p0.xcoord(),p1.xcoord()) 
#endif
                               );
  }
    
  Comparison_result compare_y(const Point_2& p0, const Point_2& p1) const {
    return (Comparison_result)(

#if (__LEDA__ >= 380)
                               // full specification is required for msvc:
                               Self::Point_2::cmp_y(p0,p1)
#else // backward compatability to LEDA
                               compare(p0.ycoord(),p1.ycoord())
#endif
                               );
  }

  //on X_curve only - not Curve!
  bool curve_is_vertical(const X_curve_2& cv) const {
    CGAL_assertion(is_x_monotone(cv));
    return compare_x(curve_source(cv),curve_target(cv))==EQUAL;
  } 

  bool curve_is_in_x_range(const X_curve_2& cv, const Point_2& p) const {
    CGAL_assertion(is_x_monotone(cv));
    return (compare_x(p,curve_source(cv)) * 
            compare_x(p,curve_target(cv))) <=0 ;
  }

  Comparison_result curve_compare_at_x(const X_curve_2& cv1, 
                                       const X_curve_2& cv2, 
                                       const Point_2& p) const 
  {
    CGAL_assertion(is_x_monotone(cv1));
    CGAL_assertion(is_x_monotone(cv2));
      
    if (!curve_is_in_x_range(cv1,p) || !curve_is_in_x_range(cv2,p) )
      return EQUAL;
      
    typename X_curve_2::const_iterator pit=cv1.begin();
    typename X_curve_2::const_iterator after=pit; ++after;
      
    for( ; (compare_x(*pit,p) * compare_x(*after,p))  >0 ; ++pit,++after) {}
      
    Segment l1(*pit,*after);

    pit=cv2.begin();
    after=pit; ++after;
    for( ; (compare_x(*pit,p) * compare_x(*after,p)) > 0 ; ++pit,++after) {}
      
    Segment l2(*pit,*after);

    // bug ??? in LEDA - 
    // cmp_segments_at_xcoord returns wrong answer if
    // cv1 (or cv2) are from right to left
    // cv1_ and cv2_ are the same as cv1 and cv2 - oriented from left to right
    if (compare_x(l1.target(), l1.source())==SMALLER)
      l1 = l1.reversal();
    if (compare_x(l2.target(), l2.source())==SMALLER)
      l2 = l2.reversal();
      
    //leda function :
    return (Comparison_result)CGAL_LEDA_SCOPE::cmp_segments_at_xcoord(l1,l2,p);

  }
  
  
  Comparison_result curve_compare_at_x_left(const X_curve_2& cv1, 
                                            const X_curve_2& cv2,
                                            const Point_2& p) const 
  {
    CGAL_assertion(is_x_monotone(cv1));
    CGAL_assertion(is_x_monotone(cv2));
      
      
    if (!curve_is_in_x_range(cv1,p) || !curve_is_in_x_range(cv2,p) )
      return EQUAL;
      
    Point_2 leftmost1 = 
      (compare_x(curve_source(cv1),curve_target(cv1))==LARGER) ? 
      curve_target(cv1) : curve_source(cv1);
    Point_2 leftmost2 = 
      (compare_x(curve_source(cv2),curve_target(cv2))==LARGER) ? 
      curve_target(cv2) : curve_source(cv2);
      
    //special cases wher returns EQUAL
    if (curve_is_vertical(cv1) || (curve_is_vertical(cv2))) return EQUAL;
    if (compare_x(leftmost1,p)!=SMALLER || compare_x(leftmost2,p)!=SMALLER) {
      return EQUAL;
    }
      
    typename X_curve_2::const_iterator pit=cv1.begin();
    typename X_curve_2::const_iterator after=pit; ++after;
    
    for( ; (compare_x(*pit,p) * compare_x(*after,p)) > 0 ; ++pit,++after) {}
      
    Segment l1(*pit,*after);

    pit=cv2.begin();
    after=pit; ++after;
    for( ; (compare_x(*pit,p) * compare_x(*after,p)) > 0 ; ++pit,++after) {}
      
    Segment l2(*pit,*after);
    
    //leda func :
    Comparison_result r = 
      (Comparison_result)CGAL_LEDA_SCOPE::cmp_segments_at_xcoord(l1,l2,p);
    
    if ( r != EQUAL)
      return r;
    else {
      //check if they are right endpoints 
      //(and compare to the one from the left) otherwise -
      if (compare_x(l1.source(),l1.target())==SMALLER)
        l1=l1.reversal();
      if (compare_x(l2.source(),l2.target())==SMALLER)
        l2=l2.reversal();

      return (Comparison_result)CGAL_LEDA_SCOPE::cmp_slopes(l2,l1);
    }

  } 


  Comparison_result curve_compare_at_x_right(const X_curve_2& cv1, 
                                             const X_curve_2& cv2,
                                             const Point_2& p) const 
  {
    CGAL_assertion(is_x_monotone(cv1));
    CGAL_assertion(is_x_monotone(cv2));
      
    if (!curve_is_in_x_range(cv1,p) || !curve_is_in_x_range(cv2,p))
      return EQUAL;
      
    Point_2 rightmost1 = 
      (compare_x(curve_source(cv1),curve_target(cv1))==SMALLER)
      ? curve_target(cv1) : curve_source(cv1);
    Point_2 rightmost2 = 
      (compare_x(curve_source(cv2),curve_target(cv2))==SMALLER)
      ? curve_target(cv2) : curve_source(cv2);
      
    //special cases where returns EQUAL
    if (curve_is_vertical(cv1) || (curve_is_vertical(cv2))) return EQUAL;
    if (compare_x(rightmost1,p)!=LARGER || compare_x(rightmost2,p)!=LARGER) {
      return EQUAL;
    }
      
    if (!curve_is_in_x_range(cv1,p) || !curve_is_in_x_range(cv2,p) )
      return EQUAL;
      
    typename X_curve_2::const_iterator pit=cv1.begin();
    typename X_curve_2::const_iterator after=pit; ++after;
      
    for( ; (compare_x(*pit,p) * compare_x(*after,p)) > 0 ; ++pit,++after) {}
      
    Segment l1(*pit,*after);
      
    pit=cv2.begin();
    after=pit; ++after;
    for( ; (compare_x(*pit,p) * compare_x(*after,p)) > 0 ; ++pit,++after) {}
      
    Segment l2(*pit,*after);
      
    Comparison_result r = 
      (Comparison_result)CGAL_LEDA_SCOPE::cmp_segments_at_xcoord(l1,l2,p);
      
    if ( r != EQUAL) {
      return r;
    }
    else
      //check if they are left endpoints (and compare to the one from the 
      //right) otherwise -
      { 
        CGAL_assertion(curve_compare_at_x(cv1,cv2,p) == EQUAL);
        //return (CGAL::compare_y_at_x(point_to_right(p),l1,l2)); 
        
        if (compare_x(l1.source(),l1.target())==LARGER)
          l1=l1.reversal();
        if (compare_x(l2.source(),l2.target())==LARGER)
          l2=l2.reversal();
        
        return (Comparison_result)CGAL_LEDA_SCOPE::cmp_slopes(l1,l2);
        
      }
  }
  
  
  
  Curve_point_status 
  curve_get_point_status(const X_curve_2 &cv, const Point_2& p) const
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
      
    typename X_curve_2::const_iterator pit=cv.begin();
    typename X_curve_2::const_iterator after=pit; ++after;
    while ( (compare_x(*pit,p) * compare_x(*after,p)) > 0 ) {
      ++pit; ++after;
    }
      
    Segment l(*pit,*after);
      
    //taken from Eyal's LedaTraits
    Orientation o;
    if (compare_x(l.source(), l.target()) < 0)
      o = (Orientation)CGAL_LEDA_SCOPE::orientation(l.source(), l.target(), p);
    else
      o = (Orientation)CGAL_LEDA_SCOPE::orientation(l.target(), l.source(), p);
      
    if (o < 0)
      return UNDER_CURVE;
    if (o > 0)
      return ABOVE_CURVE;
    if (o == 0)
      return ON_CURVE;
    return ON_CURVE;
      
  }
  
  
  //precondition - same as in pm
  bool curve_is_between_cw(const X_curve_2& cv,const X_curve_2& first,
                           const X_curve_2& second, const Point_2& p) const
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
      
      

    X_curve_2 cv0 = first;
    X_curve_2 cv1 = second;
    X_curve_2 cvx = cv;
      
      
      
    if ( !is_same(curve_source(cv0),p) ) cv0 = curve_flip(cv0);
    if ( !is_same(curve_source(cv1),p) ) cv1 = curve_flip(cv1);
    if ( !is_same(curve_source(cvx),p) ) cvx = curve_flip(cvx);
      
    typename X_curve_2::iterator xcit=cv0.begin();++xcit;
    Point_2 p0(*xcit);
    xcit=cv1.begin(); ++xcit;
    Point_2 p1(*xcit);
    xcit=cvx.begin(); ++xcit;
    Point_2 px(*xcit);
      
    if (is_same(p0,p1))
      return true; 
      
    int orient0=CGAL_LEDA_SCOPE::orientation(p0,p,px);
    int orient1=CGAL_LEDA_SCOPE::orientation(p1,p,px);
    int orient=orient0*orient1;
      
    if (orient < 0) { //one is a left_turn the other right_turn
      return (orient0 == LEFT_TURN); //left_turn
    }
    else { //both are either left or right turns (or one is colinear)
      return (CGAL_LEDA_SCOPE::orientation(p0,p,p1)==RIGHT_TURN); //right_turn
    }
  }

  bool point_is_same(const Point_2 & p, const Point_2 & q) const
  {
    return is_same(p, q);
  }
  
  bool curve_is_same(const X_curve_2& cv1, const X_curve_2& cv2) const {
    CGAL_assertion(is_x_monotone(cv1));
    CGAL_assertion(is_x_monotone(cv2));
    
    if (cv1.size()!=cv2.size())
      return false;
    typename X_curve_2::const_iterator it1=cv1.begin();
    typename X_curve_2::const_iterator it2=cv2.begin();
    while (it1!=cv1.end()) {
      if (!is_same(*it1++,*it2++))
        return false;
    }
    return true;
  }
  
  Point_2 curve_source(const X_curve_2& cv) const {
    return *(cv.begin());
  }
  
  Point_2 curve_target(const X_curve_2& cv) const {
    typename X_curve_2::const_iterator it=cv.end(); --it;
    return *it;
  }
  
  ///////////////////////////////////////////////////////
  //         ARRANGEMENT FUNCS



  X_curve_2 curve_flip(const X_curve_2& cv) const {
    X_curve_2 cv1(cv);
    std::reverse(cv1.begin(),cv1.end()); //more generic
    return cv1; 
  }

  bool is_x_monotone(const Curve_2& cv) const {
    CGAL_assertion(cv.size()>=2); //one point is not a curve

    if (cv.size()==2) return true; //segment

    typename X_curve_2::const_iterator p0=cv.begin();
    typename X_curve_2::const_iterator p1=p0; ++p1;
    typename X_curve_2::const_iterator p2=p1; ++p2;

    for(; p2!=cv.end(); ++p0,++p1,++p2) {
      if ( compare_x(*p0,*p1) * compare_x(*p1,*p2) <=0 )
        return false;
    }

    //<= a matter of decision - only one vertical segment is considered 
    //x-monotone
    return true;

  }

  //cuts into x-monotone curves, each vertical segment is 1 x-monotone curve
  //and not part of a bigger x-monotone polyline
  void make_x_monotone(const Curve_2& cv, std::list<Curve_2>& l) 
  {
    CGAL_precondition(cv.size() >= 2); //one point is not a curve
    
    // If polyline is a segment
    if (cv.size() == 2) 
    { 
      l.push_back(cv);
      return;
    }

    typename X_curve_2::const_iterator p0 = cv.begin();
    typename X_curve_2::const_iterator p1 = p0; ++p1;
    typename X_curve_2::const_iterator p2 = p1; ++p2;

    typename X_curve_2::const_iterator last_cut = p0;

    // For each consecutive triplet of points
    for(; p2 != cv.end(); ++p0, ++p1, ++p2) 
    {
      if (compare_x(*p0,*p1) == EQUAL) { //vertical segment - (p0,p1) 
        if (p0!=last_cut) { //needed for the case:
          //    | p1
          //    o p0=lastcut (was p1 before)
          //    |

          l.push_back(X_curve_2(last_cut,p1)); 
          //push_back the curve (last_cut...p0)
        }
        l.push_back(X_curve_2(p0,p2)); //push_back the segment (p0,p1)
        last_cut=p1;
      }
      else
        if ( compare_x(*p0,*p1) * compare_x(*p1,*p2) <= 0 ) {
          l.push_back(X_curve_2(last_cut,p2));
          last_cut=p1;
        }
    }

    l.push_back(X_curve_2(last_cut,p2)); //push the residue (last cut to end)


    CGAL_assertion(p2==cv.end());
    
  }

  void curve_split(const X_curve_2& cv, X_curve_2& c1, X_curve_2& c2, 
                   const Point_2& split_pt) {
    
    //split curve at split point into c1 and c2

    CGAL_precondition(curve_source(cv) != split_pt);
    CGAL_precondition(curve_target(cv) != split_pt);

    typename X_curve_2::const_iterator p0=cv.begin();
    typename X_curve_2::const_iterator p1=p0; ++p1;
 
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

    typename X_curve_2::const_iterator ci=cv.begin();
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

  }

  bool do_intersect_to_right(const X_curve_2 & ca, 
                             const X_curve_2 & cb,
                             const Point_2   & pt) const 
  {
    bool res = do_intersect_to_right_1(ca, cb, pt);

    Curve_2 cv3(ca), cv4(cb);
    cv3=curve_flip(cv3);
    cv4=curve_flip(cv4);
    CGAL_assertion(do_intersect_to_right_1(cv3, cv4, pt) == res);

    return res;
  }

  //returns true iff the intersection is lexicographically strictly right of pt

  bool do_intersect_to_right_1(const X_curve_2 & ca, 
                             const X_curve_2 & cb,
                             const Point_2   & pt) const 
  {
    CGAL_precondition(is_x_monotone(ca));
    CGAL_precondition(is_x_monotone(cb));

    // Check if both first points are left of pt, if they are reach the 
    // points directly left of pt, and check if their segments intersect  
    // to the right of pt, if not continue with a normal sweep until 
    // we find an intersection point or we reach the end.

    // Local copies of curves
    X_curve_2 c1(ca),c2(cb);

    // We make sure the local curves are oriented from left to right
    //uses LEDA's compare function that compares two rat_points 
    //lexicographically
    if (CGAL_LEDA_SCOPE::compare(curve_source(ca),curve_target(ca)) > 0 )
      c1=curve_flip(ca);
    if (CGAL_LEDA_SCOPE::compare(curve_source(cb),curve_target(cb)) > 0)
      c2=curve_flip(cb);

    typename X_curve_2::const_iterator i1s = c1.begin();
    typename X_curve_2::const_iterator i1t = i1s; ++i1t;

    typename X_curve_2::const_iterator i2s = c2.begin();
    typename X_curve_2::const_iterator i2t = i2s; ++i2t;

    // For some reason MSVC requires the following:
    typename X_curve_2::const_iterator c1End = c1.end();
    typename X_curve_2::const_iterator c2End = c2.end();

    //increment this variable if curve starts left of pt
    int number_to_left = 0; 

    if (CGAL_LEDA_SCOPE::compare(*i1s, pt) <= 0) {
      //increment to nearest from the left of pt
      ++number_to_left;
      for (; i1t != c1End; ++i1s, ++i1t) 
      {
        if (CGAL_LEDA_SCOPE::compare(*i1t, pt) > 0) break;
      }
      if (i1t == c1End) return false; //c1 ends to the left of pt
    }
    //now i1s holds the source vertex and i1t holds the target

    if (CGAL_LEDA_SCOPE::compare(*i2s, pt) <= 0) {
      ++number_to_left;
      for (; i2t != c2End; ++i2s, ++i2t) {
        if (CGAL_LEDA_SCOPE::compare(*i2t, pt) > 0) break;
      }
      if (i2t == c2End) return false; //c2 ends to the left of pt
    }

    if (number_to_left == 2) {
      // Check if intersection exists and is lexicographically larger
      Segment i_seg;
      Segment s1(*i1s, *i1t), s2(*i2s, *i2t);
      if (s1.intersection(s2, i_seg)) 
      {
        if ( (CGAL_LEDA_SCOPE::compare(i_seg.source(), pt) > 0) ||
             (CGAL_LEDA_SCOPE::compare(i_seg.target(), pt) > 0) )
          return true;
      }
    
      // advance to the nearer point
      if (CGAL_LEDA_SCOPE::compare(*i2t, *i1t) > 0) {
        ++i1s; ++i1t;
        if (i1t == c1End) return false;
      }
      else {
        ++i2s; ++i2t;
        if (i2t == c2End) return false;
      }

    }

    // Now we can start sweeping the chains
    while ( true ) 
    {
      // Check for intersection of the segments
      Segment i_seg;
      Segment s1(*i1s, *i1t), s2(*i2s, *i2t);
      if (s1.intersection(s2, i_seg)) {
        return true;
      }

      if (CGAL_LEDA_SCOPE::compare(*i2t, *i1t) > 0) 
      {
        ++i1s; ++i1t;
        if (i1t == c1End) return false;
      }
      else {
        ++i2s; ++i2t;
        if (i2t == c2End) return false;
      }
    }
  }

  bool nearest_intersection_to_right(const X_curve_2& cv1,
                                     const X_curve_2& cv2,
                                     const Point_2& pt,
                                     Point_2& p1,
                                     Point_2& p2) const 
  {
    bool res = nearest_intersection_to_right_1(cv1, cv2, pt, p1, p2);

    /*!\todo
      The following lines have to be under comments, since they add
      redundency in computation.
      However, we found a bug. Try to use the test of map-overlay
      with leda polyline traits.
      There is a function in the test computing intersection points,
      when using this function, we have that (may be the last) 
      intersection point to be computed cause a bug, 
      since the reversed curves (by flip) do not intersect.
      Hence p3 and p4 remain uninitialized.
    */
    //Point_2 p3, p4;
    //Curve_2 cv3(cv1), cv4(cv2);
    //cv3=curve_flip(cv3);
    //cv4=curve_flip(cv4);
    //nearest_intersection_to_right_1(cv3, cv4, pt, p3, p4);

    //cout<<"p1="<<p1<<"  p3="<<p3<<endl;
    //cout<<"p2="<<p1<<"  p4="<<p3<<endl;

    //CGAL_assertion( p1 == p3 );
    //CGAL_assertion( p2 == p4 );
    
    return res;
  }

  // NOTE: in the following function:
  // if there is an overlap, we return only a segment of overlap 
  // (not the whole chain) - this is easier to implement and is sufficient 
  // for our needs
  bool nearest_intersection_to_right_1(const X_curve_2& cv1,
                                     const X_curve_2& cv2,
                                     const Point_2& pt,
                                     Point_2& p1,
                                     Point_2& p2) const 
  {
    CGAL_assertion(is_x_monotone(cv1));
    CGAL_assertion(is_x_monotone(cv2));

    bool found = false;

    if ( ! do_intersect_to_right(cv1,cv2,pt)) return false;
    
    // Makes x increase.
    X_curve_2 c1(cv1),c2(cv2);
    if (CGAL_LEDA_SCOPE::compare(curve_source(c1),curve_target(c1)) > 0)
      c1=curve_flip(cv1);
    if (CGAL_LEDA_SCOPE::compare(curve_source(c2),curve_target(c2)) > 0)
      c2=curve_flip(cv2);

    // check if both first points are left of pt, if they are reach the 
    // points directly left of pt, and check if their segments intersect  
    //to the right of pt, if not continue with a normal sweep until 
    //we find an intersection point or we reach the end.
    
    typename X_curve_2::const_iterator i1s=c1.begin(), i1e = c1.end();
    typename X_curve_2::const_iterator i1t=i1s; ++i1t;

    typename X_curve_2::const_iterator i2s=c2.begin(), i2e = c2.end();
    typename X_curve_2::const_iterator i2t=i2s; ++i2t;

    // For some reason MSVC requires the following:
    typename X_curve_2::const_iterator c1End = c1.end();
    typename X_curve_2::const_iterator c2End = c2.end();

    int number_to_left=0; //increment this variable if curve starts left of pt
    
    if (CGAL_LEDA_SCOPE::compare(*i1s, pt) <= 0) {
      //increment to nearest from the left of pt
      ++number_to_left;
      for (; i1t != i1e; ++i1s,++i1t) {
        if (CGAL_LEDA_SCOPE::compare(*i1t, pt) > 0) break;
      }
      CGAL_assertion (i1t != c1End);
    }

    if (CGAL_LEDA_SCOPE::compare(*i2s, pt) <= 0) {
      ++number_to_left;
      for (; i2t != i2e; ++i2s, ++i2t) {
        if (CGAL_LEDA_SCOPE::compare(*i2t, pt) > 0) break;
      }
      CGAL_assertion(i2t != c2End)  ; //c2 ends to the left of pt
    }
   
    if (number_to_left == 2) {
      //check if intersection exists and is lexicographically larger
      Segment s1(*i1s,*i1t),s2(*i2s,*i2t);
      Segment i_seg;
      if( s1.intersection(s2, i_seg)) {
        if( CGAL_LEDA_SCOPE::compare( i_seg.target(),i_seg.source() ) == 0 )
        {
          // Intersection is only one point
          if( CGAL_LEDA_SCOPE::compare( i_seg.target(), pt ) > 0 )
          {
            p1 = p2 = point_normalize(i_seg.target());
            found = true;
          }
        }
        if( CGAL_LEDA_SCOPE::compare( i_seg.target(),i_seg.source() ) != 0 &&
            ! found ) 
        {
          // Intersection is a segment
          if( CGAL_LEDA_SCOPE::compare( i_seg.target(), i_seg.source()) < 0 )
            i_seg = i_seg.reverse();
          if( CGAL_LEDA_SCOPE::compare( i_seg.target(),pt ) > 0 ) {
            p2=i_seg.target();
            if( CGAL_LEDA_SCOPE::compare( i_seg.source(),pt ) > 0 ){
              p1=i_seg.source();                           
            }
            else {
              //p1=pt;
              // Modified by Eug
              // Performing vertical ray shooting from pt.
              // Finding the intersection point. We know by now
              // that there is exactly ONE point. Assinging this 
              // point to p1.
              Point_2 ap1( pt.xcoord(), i_seg.source().ycoord() );
              Point_2 ap2( pt.xcoord(), i_seg.target().ycoord() );
              Segment vertical_pt_x_base( ap1, ap2 );
              Segment i_vertical;
              vertical_pt_x_base.intersection( i_seg, i_vertical );
              p1 = i_vertical.source();
            }
            found = true;
          }
        }
      }//endif intersection 
      if ( ! found) {
        //advance to the nearer point
        if( CGAL_LEDA_SCOPE::compare( *i2t, *i1t ) > 0 ) {
          ++i1s; ++i1t;
          if( i1t == i1e ) 
            return false;
        }
        else {
          ++i2s; ++i2t;
          if( i2t == i2e ) 
            return false;
        }
      }
    }//endif number_to_left==2
 
    while( !found ){
      //check intersection of the segments
      Segment s1(*i1s,*i1t),s2(*i2s,*i2t);
      Segment i_seg;
       
      if( s1.intersection( s2,i_seg ) ){
        if( CGAL_LEDA_SCOPE::compare( i_seg.target(),i_seg.source() ) == 0 ){
          // Intersection is only one point
          if( CGAL_LEDA_SCOPE::compare( i_seg.target(), pt ) > 0 ){
            p1 = p2 = point_normalize(i_seg.target());
            found = true;
          }
        }
        if( CGAL_LEDA_SCOPE::compare( i_seg.target(),i_seg.source() ) != 0  &&
            !found ){
          //check if intersection seg to the right of pt
          if( CGAL_LEDA_SCOPE::compare( i_seg.source(), i_seg.target() ) > 0 )
            i_seg=i_seg.reverse();

          if( CGAL_LEDA_SCOPE::compare( i_seg.target(), pt ) > 0 ){
            p2 = i_seg.target();
            if( CGAL_LEDA_SCOPE::compare( i_seg.source(), pt ) > 0 ){
              p1=i_seg.source();
            }
            else {
              // p1=pt; 
              // Modified by Eug
              // Performing vertical ray shooting from pt.
              // Finding the intersection point. We know by now
              // that there is exactly ONE point. Assinging this 
              // point to p1. 
              Point_2 ap1( pt.xcoord(), i_seg.source().ycoord() );
              Point_2 ap2( pt.xcoord(), i_seg.target().ycoord() );
              Segment vertical_pt_x_base( ap1, ap2 );
              Segment i_vertical;
              vertical_pt_x_base.intersection( i_seg, i_vertical );
              p1 = i_vertical.source();
            }           
            found = true;
          }
        }
      } // endif s1.intersection  
      if( !found ){
        //advance to the nearer point
        if( CGAL_LEDA_SCOPE::compare( *i2t, *i1t ) > 0 ){
          ++i1s; ++i1t;
          if( i1t == i1e ) 
            return false;
        }
        else {
          ++i2s; ++i2t;
          if( i2t == i2e ) 
            return false;
        }
      } 
    } // end while ( ! found)
       
    // if the x point is at the end of a segment, then there might be 
    // an overlap in the continious of the polyline
    if ( found && CGAL_LEDA_SCOPE::compare( p1, p2 ) == 0 ){
      typename X_curve_2::const_iterator s1=i1s, t1=i1t, s2=i2s, t2=i2t;
      s1++; t1++; s2++; t2++;
      if( t1 != i1e && t2 != i2e){

        // check for overlap after x point
        Segment i_seg;

        if ( CGAL_LEDA_SCOPE::compare( p1, *i1t ) == 0  && 
             CGAL_LEDA_SCOPE::compare( p1, *i2t ) == 0 ){
          Segment tmp_seg1( *s1, *t1 );
          Segment tmp_seg2( *s2, *t2 );
	  if( !tmp_seg1.intersection( tmp_seg2, i_seg ) )
            return false;
        }
        else if ( CGAL_LEDA_SCOPE::compare( p1, *i1t ) == 0 ){
          Segment tmp_seg1( *s1, *t1 );
          Segment tmp_seg2( *i2s, *i2t );
	  if( !tmp_seg1.intersection( tmp_seg2, i_seg ) )
            return false;
        }
        else if ( is_same( p1, *i2t)) {
          Segment tmp_seg1( *i1s, *i1t );
          Segment tmp_seg2( *s2, *t2 );
	  if( !tmp_seg1.intersection( tmp_seg2, i_seg ) )
            return false;
        }
        // no need to check whether intersection seg to the right of pt, 
        // because we have already found an x point to the right of pt.
        if( CGAL_LEDA_SCOPE::compare( i_seg.source(), i_seg.target() ) > 0 )
          i_seg=i_seg.reverse();
        p1 = i_seg.source();
        p2 = i_seg.target();
      }        
    }
    return found;
  }

  bool curves_overlap(const X_curve_2& ca, const X_curve_2& cb) const 
  {
    CGAL_assertion(curves_overlap_1(ca, cb) == 
                   curves_overlap_1(curve_flip(ca), curve_flip(cb)));
    return curves_overlap_1(ca, cb);
  }
  bool curves_overlap_1(const X_curve_2& ca, const X_curve_2& cb) const 
  {
    CGAL_assertion(is_x_monotone(ca));
    CGAL_assertion(is_x_monotone(cb));

    // c1 and c2 are local copies of curves
    X_curve_2 c1(ca), c2(cb);

    // Make sure local copies are oriented from left to right 
    // We use LEDA's compare() function that compares two rat_points 
    // lexicographically.
    if (CGAL_LEDA_SCOPE::compare(curve_source(ca), curve_target(ca)) > 0 )
      c1=curve_flip(ca);
    if (CGAL_LEDA_SCOPE::compare(curve_source(cb), curve_target(cb)) > 0 )
      c2=curve_flip(cb);

    // i1s and i1t are iterators over the points of the polylines.
    // i1t is always one point further after i1s, so if we start with
    // i1s == c1.begin() and increment i1s and i1t together, (*i1s, *i1t)
    // is the current segment of the polyline.
    typename X_curve_2::const_iterator i1s = c1.begin();
    typename X_curve_2::const_iterator i1t = i1s; ++i1t;

    typename X_curve_2::const_iterator i2s = c2.begin();
    typename X_curve_2::const_iterator i2t = i2s; ++i2t;

    // For some reason MSVC requires the following:
    typename X_curve_2::const_iterator c1End = c1.end();
    typename X_curve_2::const_iterator c2End = c2.end();

    // Now i1s holds the source vertex and i1t holds the target
    Segment i_seg;
    Segment s1(*i1s, *i1t), s2(*i2s, *i2t);
    if (s1.intersection(s2,i_seg)) 
    {
      // If the resulting intersection segment (can be 0 lengthed) is not
      // of zero length then there is an overlap.
      if (CGAL_LEDA_SCOPE::compare(i_seg.source(), i_seg.target()) != 0)
        return true;
    }
    
    // Advance to the nearer point
    if (CGAL_LEDA_SCOPE::compare(*i2t,*i1t) > 0) {
      ++i1s; ++i1t;
      if (i1t == c1End) return false;
    }
    else {
      ++i2s; ++i2t;
      if (i2t == c2End) return false;
    }

    //NOW we can start sweeping the chains

    while ( true ) 
    {
      Segment i_seg;
      Segment s1(*i1s, *i1t), s2(*i2s, *i2t);
      if (s1.intersection(s2,i_seg)) 
      {
        if (CGAL_LEDA_SCOPE::compare(i_seg.source(),i_seg.target()) != 0)
          return true;
      }

      if (CGAL_LEDA_SCOPE::compare(*i2t,*i1t) > 0) {
        ++i1s; ++i1t;
        if (i1t == c1End) return false;
      }
      else {
        ++i2s; ++i2t;
        if (i2t == c2End) return false;
      }
    }
  }

  X_curve_2 curve_reflect_in_x_and_y( const X_curve_2& cv) const
  {
    X_curve_2 reflected_cv;
    typename Curve_2::const_iterator it  = cv.begin();
    for  (; it != cv.end(); it++)
      {
        reflected_cv.push_back( point_reflect_in_x_and_y( *it)); 
      }
    return reflected_cv;
  }


  Point_2 point_reflect_in_x_and_y (const Point_2& pt) const
  {
    // use hx(), hy(), hw() in order to support both Homogeneous and Cartesian
    Point_2 reflected_pt( -pt.X(), -pt.Y(), pt.W());
    return reflected_pt;
  }

  /////////////////////////////////////////////////////////////////
private:

  bool is_same(const Point_2 &p1, const Point_2 &p2) const
  {
    return (compare_x(p1, p2) == EQUAL) &&
      (compare_y(p1, p2) == EQUAL);
  }

  Point_2 point_to_left(const Point_2& p) const {
    return Point_2(p.xcoord()-1, p.ycoord());  
  }
  Point_2 point_to_right(const Point_2& p) const {
    return Point_2(p.xcoord()+1, p.ycoord());
  }
  
protected:
  Point_2 point_normalize(const Point_2 &pt) const {

    leda_integer g, x, y, w;
    x = pt.X();
    y = pt.Y();
    w = pt.W();
    if (x.iszero() &&  y.iszero()) {
      //g = w;
      return Point_2(x,y,leda_integer(1));
    }
    else {
      g = gcd(x, y);
      g = gcd(g, w);
      //x,y,w now hold zero
      
      return Point_2(x/g,y/g,w/g);
    }

  }
};

CGAL_END_NAMESPACE

#endif // CGAL_ARR_LEDA_POLYLINE_TRAITS_H
