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
// release       : 
// release_date  : 2000, September 14
//
// file          : include/CGAL/point_set_traits_2.h
// package       : Point_set_2 (1.2.4)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.2.4
// revision_date : 14 September 2000 
// author(s)     : Kurt Mehlhorn, Stefan Naeher, Matthias Baesken
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ======================================================================

#ifndef POINT_SET_2_TRAITS
#define POINT_SET_2_TRAITS

#if !defined(LEDA_ROOT_INCL_ID)
#define LEDA_ROOT_INCL_ID 400898
#include <LEDA/REDEFINE_NAMES.h>
#endif

#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Direction_2.h>

CGAL_BEGIN_NAMESPACE


#if defined(LEDA_PREFIX)
#undef vector
#endif



template<class T>
class CG_Lex_xy : public leda_cmp_base<T> {
  int operator()(const T& p1,const T& p2) const
  { 
   Comparison_result c = compare_lexicographically_xy(p1,p2);    
   switch(c){
   case SMALLER: return -1;
   case LARGER: return 1;
   case EQUAL: return 0;
   }
   return 0;
  } 
};

template<class T>
struct CG_compare_dist {
  Comparison_result operator()(const T& p1,const T& p2,const T& p3) const
  { 
    Comparison_result c= cmp_dist_to_point(p1,p2,p3); 
    return c;
  }
};


template<class T>
struct CG_orientation {
  Orientation operator()(const T& p1,const T& p2,const T& p3) const
  { 
    Orientation or=orientation(p1,p2,p3); 
    return or;
  } 
};

template<class T>
struct CG_side_of_circle {
  int operator()(const T& p1,const T& p2,const T& p3, const T& p4) const
  { 
   Oriented_side s= side_of_oriented_circle(p1,p2,p3,p4); 
   return s;
  }
};

template<class T>
struct CG_side_of_halfspace {
  typedef typename T::R  Repres;

  Orientation operator()(const T& p1,const T& p2, const T& p3) const
  { 
    Segment_2<Repres>   s(p1,p2);
    Direction_2<Repres> d = s.direction();
    
    Vector_2<Repres>    v(d.dx(),d.dy());
    Vector_2<Repres>    v2= v.perpendicular(CLOCKWISE);
    Point_2<Repres>     ph = p1 + v2;
    
    Orientation or=orientation(p1,ph,p3); 
    return or;
  }
};


template<class S, class T>
struct CG_segment_contains {
  bool operator()(const S& seg, const T& p) const
  { return seg.has_on(p); }
};

template<class T, class Rep>
struct CG_sqr_dist {
  Rep operator()(const T& p1,const T& p2) const
  { return squared_distance(p1,p2); }
};


template<class L,class T, class Rep>
struct CG_line_sqr_dist {
  Rep operator()(const L& l,const T& p) const
  { return squared_distance(l,p); }
};


template<class C,class T>
struct CG_circle_ptori {
  Bounded_side operator()(const C& c, const T& p) const
  { 
   Bounded_side s = c.bounded_side(p); 
   return s;
  }
};


template<class C,class T>
struct CG_circle_center {
  T operator()(const C& c) const
  { return c.center(); }
};


template<class C,class T>
struct CG_create_circle_ppp {
  C operator()(const T& p1,const T& p2,const T& p3) const
  { return C(p1,p2,p3); }
};

template<class C,class T>
struct CG_create_segment_pp {
  C operator()(const T& p1,const T& p2) const
  { return C(p1,p2); }
};

template<class C,class T>
struct CG_create_line_pp {
  C operator()(const T& p1,const T& p2) const
  { return C(p1,p2); }
};


template<class T>
class point_set_traits_2 {
 public:
  typedef typename T::FT                      FT;
  typedef Point_2<T>                          Point;
  typedef Circle_2<T>                         Circle;
  typedef Segment_2<T>                        Segment;
  typedef Line_2<T>                           Line;

  typedef CG_Lex_xy<Point>                    Compare_xy_2;
  typedef CG_compare_dist<Point>              Compare_dist_2;
  typedef CG_orientation<Point>               Orientation;
  typedef CG_side_of_circle<Point>            Side_of_oriented_circle_2;
  typedef CG_side_of_halfspace<Point>         Side_of_halfspace_2;
  typedef CG_segment_contains<Segment,Point>  Segment_has_on_2;
  typedef CG_sqr_dist<Point,FT>               Squared_distance;
  typedef CG_line_sqr_dist<Line,Point,FT>     Squared_distance_to_line;
  typedef CG_circle_ptori<Circle,Point>       Circle_bounded_side_2;
  typedef CG_circle_center<Circle,Point>      Circle_center_2;
  
  typedef CG_create_circle_ppp<Circle,Point>  Construct_circle_2; 
  typedef CG_create_segment_pp<Segment,Point> Construct_segment_2;
  typedef CG_create_line_pp<Line,Point>       Construct_line_2;   

  Orientation  get_orientation_object() const
  { return Orientation(); }

  Side_of_oriented_circle_2 get_side_of_oriented_circle_2_object() const
  { return Side_of_oriented_circle_2(); }

  Side_of_halfspace_2 get_side_of_halfspace_2_object() const
  { return Side_of_halfspace_2(); }

  Segment_has_on_2  get_segment_has_on_2_object() const
  { return Segment_has_on_2(); }

  Squared_distance get_squared_distance_object() const
  { return Squared_distance(); }

  Squared_distance_to_line  get_squared_distance_to_line_object() const
  { return Squared_distance_to_line(); }

  Circle_bounded_side_2 get_circle_bounded_side_2_object() const
  { return Circle_bounded_side_2(); }

  Circle_center_2 get_circle_center_2_object() const
  { return Circle_center_2(); }
  
  Construct_circle_2 get_construct_circle_2_object() const
  { return Construct_circle_2(); }
  
  Construct_segment_2 get_construct_segment_2_object() const
  { return Construct_segment_2(); }
  
  Construct_line_2 get_construct_line_2_object() const
  { return Construct_line_2(); }
  

  point_set_traits_2() { }
};

CGAL_END_NAMESPACE


#if LEDA_ROOT_INCL_ID == 400898
#undef LEDA_ROOT_INCL_ID
#include <LEDA/UNDEFINE_NAMES.h>
#endif

#endif
