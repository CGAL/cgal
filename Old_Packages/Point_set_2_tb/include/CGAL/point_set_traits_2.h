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
// release_date  : 
//
// file          : include/CGAL/point_set_traits_2.h
// package       : Point_set_2_tb (0.9)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9
// revision_date : 27 March 2001
// author(s)     : Kurt Mehlhorn, Stefan Naeher, Matthias Baesken
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ======================================================================

#ifndef POINT_SET_2_TRAITS
#define POINT_SET_2_TRAITS

#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Direction_2.h>

CGAL_BEGIN_NAMESPACE

template<class T>
struct CG_compare_dist {
  Comparison_result operator()(const T& p1,const T& p2,const T& p3) const
  { 
    Comparison_result c= cmp_dist_to_point(p1,p2,p3); 
    return c;
  }
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


// when appropriate kernel traits functionality becomes available,
// this class might be removed ...

template<class T>
class point_set_traits_2 {
 public:
  typedef Point_2<T>                          Point;
  typedef Circle_2<T>                         Circle;

  typedef CG_compare_dist<Point>              Compare_dist_2;
  typedef CG_circle_ptori<Circle,Point>       Circle_bounded_side_2;
  typedef CG_circle_center<Circle,Point>      Circle_center_2;

  point_set_traits_2() { }
};

CGAL_END_NAMESPACE



#endif
