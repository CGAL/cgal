// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// source        : 
// file          : include/CGAL/_test_cls_triangulation_3.C
// revision      : 
// revision_date : 

// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)

// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

// Points file for _test_cls_triangulation.C

#include "_test_types.h"
#include <list>

typedef Triangulation_test_point Point;
typedef list<Point>              list_point;


// list of Points for T1_0 :


Point p[5]=
    { Point(0,0,0),Point(0,0,1),Point(0,0,2),Point(0,0,3),Point(0,0,4)};
  list_point l;
  int i;
  for (i=0;i<5;i++)
    l.push_back(p[i]);

// Points for T2_0 :
  Point p1=Point(5,5,0); 
  Point p2=Point(4,4,0);
  Point p3=Point(6,6,0); // 1- dimensional until this point
  Point p4=Point(5,3,0); // 2-dimensional
  Point p5=Point(5,7,0); 
  Point p6=Point(5,4,0); 
  Point p7=Point(5,6,0); 
  Point p8=Point(0,0,0); 
  Point p9=Point(5,5,0); 

// Points for T3_1 :
Point q[22] = 
{
 Point(0,0,0) Point(4,4,0), Point(0,4,0), Point(4,0,0),
 Point(1,3,1), Point(3,1,1), Point(3,3,1), Point(1,1,1), Point(2,2,2),
 Point(1,3,3), Point(3,1,3), Point(3,3,3), Point(1,1,3), 
 Point(0,0,4) Point(4,4,4), Point(0,4,4), Point(4,0,4),
 Point(1,3,5), Point(3,1,5), Point(3,3,5), Point(1,1,5), Point(2,2,6)};

// Points for T3_2 :

list_point lp;
int a, b, c;
  for (a=1;a!=11;a++)
    for (b=1;b!=101;b++)
      for (c=1;c!=101;c++)
       lp.push_back(Point(a*b-c*a,a-b+c,a*a-c*c+b));


