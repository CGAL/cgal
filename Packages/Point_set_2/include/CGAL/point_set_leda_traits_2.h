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
// file          : include/CGAL/point_set_leda_traits_2.h
// package       : Point_set_2 (1.2.4)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.2.4
// revision_date : 14 September 2000 
// author(s)     : Kurt Mehlhorn, Stefan Naeher, Matthias Baesken
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ======================================================================

#ifndef POINT_SET_2_LEDA_TRAITS
#define POINT_SET_2_LEDA_TRAITS

#if !defined(LEDA_ROOT_INCL_ID)
#define LEDA_ROOT_INCL_ID 400899
#include <LEDA/REDEFINE_NAMES.h>
#endif

#include <LEDA/point.h>
#include <LEDA/circle.h>
#include <LEDA/segment.h>
#include <LEDA/line.h>

#include <LEDA/rat_point.h>
#include <LEDA/rat_circle.h>
#include <LEDA/rat_segment.h>
#include <LEDA/rat_line.h>

CGAL_BEGIN_NAMESPACE

class Leda_compare_ratpoints : public leda_cmp_base<leda_rat_point> {
  int operator()(const leda_rat_point& p1,const leda_rat_point& p2) const
  { 
    return leda_rat_point::cmp_xy(p1,p2);
  }
};

class Leda_compare_points : public leda_cmp_base<leda_point> {
  int operator()(const leda_point& p1,const leda_point& p2) const
  { 
    return leda_point::cmp_xy(p1,p2);
  }
};

template<class T>
struct Leda_compare_dist {
  Comparison_result operator()(const T& p1,const T& p2,const T& p3) const
  { 
    int c = p1.cmp_dist(p2,p3); 
    switch (c){
     case -1: return SMALLER;
     case 0 : return EQUAL;
    }
    return LARGER;
  }
};

template<class T>
struct Leda_orientation {
  Orientation operator()(const T& p1,const T& p2,const T& p3) const
  { 
    int o = ::orientation(p1,p2,p3); 
    switch(o){
     case -1: return RIGHTTURN;
     case  0: return COLLINEAR;
    }
    return LEFTTURN;
  } 
};

template<class T>
struct Leda_side_of_circle {
  Oriented_side operator()(const T& p1,const T& p2,const T& p3, const T& p4) const
  { 
    int s = side_of_circle(p1,p2,p3,p4); 
    switch(s){
     case -1: return ON_NEGATIVE_SIDE;
     case  0: return ON_ORIENTED_BOUNDARY;
    }
    return ON_POSITIVE_SIDE;
  }
};

template<class T>
struct Leda_side_of_halfspace {
  Orientation operator()(const T& p1,const T& p2, const T& p3) const
  { 
    int o = side_of_halfspace(p1,p2,p3); 
    switch(o){
     case -1: return RIGHTTURN;
     case  0: return COLLINEAR;    
    }
    return LEFTTURN;
  }
};

template<class S, class T>
struct Leda_segment_contains {
  bool operator()(const S& seg, const T& p) const
  { return seg.contains(p); }
};

template<class T, class Rep>
struct Leda_sqr_dist {
  Rep operator()(const T& p1,const T& p2) const
  { return p1.sqr_dist(p2); }
};

template<class L,class T, class Rep>
struct Leda_line_sqr_dist {
  Rep operator()(const L& l,const T& p) const
  { return l.sqr_dist(p); }
};

template<class C,class T>
struct Leda_circle_ptori {
  Bounded_side operator()(const C& c, const T& p) const
  { int wt = c.orientation();
    int so = c.side_of(p);
    int bs = wt*so; // <0  outside, ==0 on the circle, >0 inside
    switch (bs){
     case -1: return ON_UNBOUNDED_SIDE;
     case  0: return ON_BOUNDARY;
    }
    return ON_BOUNDED_SIDE;
  }
};

template<class C,class T>
struct Leda_circle_center {
  T operator()(const C& c) const
  { return c.center(); }
};


template<class C,class T>
struct Leda_create_circle_ppp {
  C operator()(const T& p1,const T& p2,const T& p3) const
  { return C(p1,p2,p3); }
};

template<class C,class T>
struct Leda_create_segment_pp {
  C operator()(const T& p1,const T& p2) const
  { return C(p1,p2); }
};

template<class C,class T>
struct Leda_create_line_pp {
  C operator()(const T& p1,const T& p2) const
  { return C(p1,p2); }
};


class point_set_leda_float_traits_2 {
 public:
  typedef double                              FT;
  typedef leda_point                            Point;
  typedef leda_circle                           Circle;
  typedef leda_segment                          Segment;
  typedef leda_line                             Line;
  
  typedef Leda_compare_points                   Compare_xy_2;
  typedef Leda_compare_dist<Point>              Compare_dist_2;
  typedef Leda_orientation<Point>               Orientation;
  typedef Leda_side_of_circle<Point>            Side_of_oriented_circle_2;
  typedef Leda_side_of_halfspace<Point>         Side_of_halfspace_2;
  typedef Leda_segment_contains<Segment,Point>  Segment_has_on_2;
  typedef Leda_sqr_dist<Point,FT>               Squared_distance;
  typedef Leda_line_sqr_dist<Line,Point,FT>     Squared_distance_to_line;
  typedef Leda_circle_ptori<Circle,Point>       Circle_bounded_side_2;
  typedef Leda_circle_center<Circle,Point>      Circle_center_2;
  
  typedef Leda_create_circle_ppp<Circle,Point>  Construct_circle_2; 
  typedef Leda_create_segment_pp<Segment,Point> Construct_segment_2; 
  typedef Leda_create_line_pp<Line,Point>       Construct_line_2;   

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

  point_set_leda_float_traits_2() { }
};

class point_set_leda_rat_traits_2 {
 public:
  typedef leda_rational                           FT;
  typedef leda_rat_point                          Point;
  typedef leda_rat_circle                         Circle;
  typedef leda_rat_segment                        Segment;
  typedef leda_rat_line                           Line;

  typedef Leda_compare_ratpoints                  Compare_xy_2;
  typedef Leda_compare_dist<Point>                Compare_dist_2;
  typedef Leda_orientation<Point>                 Orientation;
  typedef Leda_side_of_circle<Point>              Side_of_oriented_circle_2;
  typedef Leda_side_of_halfspace<Point>           Side_of_halfspace_2;
  typedef Leda_segment_contains<Segment,Point>    Segment_has_on_2;
  typedef Leda_sqr_dist<Point,FT>                 Squared_distance;
  typedef Leda_line_sqr_dist<Line,Point,FT>       Squared_distance_to_line;
  typedef Leda_circle_ptori<Circle,Point>         Circle_bounded_side_2;
  typedef Leda_circle_center<Circle,Point>        Circle_center_2;
  
  typedef Leda_create_circle_ppp<Circle,Point>    Construct_circle_2; 
  typedef Leda_create_segment_pp<Segment,Point>   Construct_segment_2;
  typedef Leda_create_line_pp<Line,Point>         Construct_line_2;   


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
  

  point_set_leda_rat_traits_2() { }
};

CGAL_END_NAMESPACE


#if LEDA_ROOT_INCL_ID == 400899
#undef LEDA_ROOT_INCL_ID
#include <LEDA/UNDEFINE_NAMES.h>
#endif

#endif
