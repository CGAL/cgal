// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       :
// release_date  : 
//
// file          : include/CGAL/Convex_hull_projective_xy_traits_2.h
// package       : Convex_hull_2 
// revision      : $Revision$
// revision_date : $Date#
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CONVEX_HULL_PROJECTIVE_XY_TRAITS_2_H
#define CONVEX_HULL_PROJECTIVE_XY_TRAITS_2_H

#include <CGAL/predicates/kernel_ftC2.h>

namespace CGAL {

template <class Point_3>
class Less_xy_plane_xy_2 
{
public:
   bool 
   operator()(const Point_3& p, const Point_3& q) const
   { 
     return 
        compare_lexicographically_xyC2(p.x(), p.y(), q.x(), q.y()) == SMALLER;
   }
};

template <class Point_3>
class Less_yx_plane_xy_2 
{
public:
   bool 
   operator()(const Point_3& p, const Point_3& q) const
   { 
     return 
        compare_lexicographically_xyC2(p.y(), p.x(), q.y(), q.x()) == SMALLER;
   }
};

template <class Point_3>
class Leftturn_plane_xy_2 
{
public:
   bool 
   operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
   { 
    return orientationC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y()) == LEFTTURN;
   }
};

template <class Point_3>
class Left_of_line_plane_xy_2
{
public:
   Left_of_line_plane_xy_2(const Point_3& a, const Point_3& b):
      p_a(a), p_b(b) 
   { }

   bool 
   operator()(const Point_3& c) const
   {
      return orientationC2(p_a.x(), p_a.y(), p_b.x(), p_b.y(), c.x(), c.y()) ==
             LEFTTURN; 
   }
private:
   Point_3 p_a;
   Point_3 p_b;
};

template <class Point_3>
class Less_dist_to_line_plane_xy_2
{
public:
   Less_dist_to_line_plane_xy_2(const Point_3& a, const Point_3& b):
      p(a), q(b)  
   { }

   bool
   operator()(const Point_3& r, const Point_3& s) const
   {
      Comparison_result
         res = cmp_signed_dist_to_lineC2(p.x(), p.y(), q.x(), q.y(),
                                         r.x(), r.y(), s.x(), s.y());
      if ( res == LARGER )
         return false;
      else if ( res == SMALLER )
         return true;
      else
         return compare_lexicographically_xyC2(r.x(), r.y(), s.x(), s.y()) 
             == SMALLER;
   }
private:
   Point_3 p;
   Point_3 q;
};

template <class Point_3>
class Less_rotate_ccw_plane_xy_2
{
public:
   Less_rotate_ccw_plane_xy_2(const Point_3& p):
      rot_point(p)
   { }

   bool
   operator()(const Point_3& p, const Point_3& q) const
   {
      Orientation orient =
               orientationC2(rot_point.x(), rot_point.y(), p.x(), p.y(), 
                             q.x(), q.y());
      if ( orient ==  LEFTTURN )
         return true;
      else if ( orient == RIGHTTURN )
         return false;
      else
      {
         if (p.x() == rot_point.x() && p.y() == rot_point.y()) return false;
         if (q.x() == rot_point.x() && q.y() == rot_point.y()) return true;
         if (p.x() == q.x() && p.y() == q.y()) return false;
         return 
            collinear_are_ordered_along_lineC2(rot_point.x(), rot_point.y(), 
                                               q.x(), q.y(), p.x(), p.y()); 
      }
   }

private:
   Point_3 rot_point;
};


template <class Point_3>
class Convex_hull_projective_xy_traits_2 
{
public:
    typedef Point_3                             Point_2;
    typedef Less_xy_plane_xy_2<Point_3>         Less_xy_2;
    typedef Less_yx_plane_xy_2<Point_3>         Less_yx_2;
    typedef Leftturn_plane_xy_2<Point_3>        Leftturn_2;
    typedef Left_of_line_plane_xy_2<Point_3>    Left_of_line_2;
    typedef Less_rotate_ccw_plane_xy_2<Point_3> Less_rotate_ccw_2;
    typedef Less_dist_to_line_plane_xy_2<Point_3> 
                                                Less_signed_distance_to_line_2;

    Left_of_line_2
    left_of_line_2_object(Point_2 p, Point_2 q) const
    {  return Left_of_line_2(p, q); }

    Less_xy_2
    less_xy_2_object() const
    {  return Less_xy_2(); }

    Less_yx_2
    less_yx_2_object() const
    {  return Less_yx_2(); }

    Leftturn_2
    leftturn_2_object() const
    {  return Leftturn_2(); }

    Less_rotate_ccw_2
    less_rotate_ccw_2_object() const
    {  return Less_rotate_ccw_2(); }

    Less_signed_distance_to_line_2
    less_signed_distance_to_line_2_object() const
    {  return Less_signed_distance_to_line_2(); }
};

}
#endif // CONVEX_HULL_PROJECTIVE_XY_TRAITS_2_H
