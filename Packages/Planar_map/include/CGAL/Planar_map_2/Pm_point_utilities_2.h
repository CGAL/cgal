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
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/Planar_map_2/Pm_point_utilities_2.h
// package       : Planar_map (5.87)
// maintainer    : Efi Fogel <efif@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Shai Hirsch       <shaihi@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_POINT_UTILITIES_2_H
#define CGAL_PM_POINT_UTILITIES_2_H

CGAL_BEGIN_NAMESPACE

template<class Point_2>
bool is_left(const Point_2 & p1, const Point_2 & p2)
{
  typedef typename Point_2::R Kernel;
  return Kernel().less_x_2_object()(p1, p2); 
}


template<class Point_2>
bool is_right(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return Kernel().less_x_2_object()(p2, p1);
}

template<class Point_2>
bool is_same_x(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return Kernel().equal_x_object()(p1, p2);
}

template<class Point_2>
bool is_lower(const Point_2 & p1, const Point_2 & p2)
{
  typedef typename Point_2::R Kernel;
  return Kernel().less_y_2_object()(p1, p2);
}

template<class Point_2>
bool is_higher(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return Kernel().less_y_2_object()(p2, p1);
}

template<class Point_2>
bool is_same_y(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return Kernel().equal_y_object()(p1, p2);
}

template<class Point_2>
bool is_same(const Point_2 & p1, const Point_2 &p2)
{
  typedef typename Point_2::R Kernel;
  return (compare_x(p1, p2) == EQUAL) && (compare_y(p1, p2) == EQUAL);
}

template<class Point_2>
const Point_2 & leftmost(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return (is_left(p1, p2) ? p1 : p2); 
}

template<class Point_2>
const Point_2 & rightmost(const Point_2 & p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return (is_right(p1, p2) ? p1 : p2); 
}
  
template<class Point_2>
const Point_2 & lowest(const Point_2 &p1, const Point_2 & p2)
{ 
  typedef typename Point_2::R Kernel;
  return (is_lower(p1, p2) ? p1 : p2);
}
  
template<class Point_2>
const Point_2 & highest(const Point_2 &p1, const Point_2 & p2)
{
  typedef typename Point_2::R Kernel;
  return (is_higher(p1, p2) ? p1 : p2);
}

CGAL_END_NAMESPACE

#endif // CGAL_PM_POINT_UTILITIES_2_H
// EOF
