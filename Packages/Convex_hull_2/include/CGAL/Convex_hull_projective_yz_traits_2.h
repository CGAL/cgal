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
// file          : include/CGAL/Convex_hull_projective_yz_traits_2.h
// package       : Convex_hull_2 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Susan Hert
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CONVEX_HULL_PROJECTIVE_YZ_TRAITS_2_H
#define CONVEX_HULL_PROJECTIVE_YZ_TRAITS_2_H

#include <CGAL/predicates/kernel_ftC2.h>
#include <CGAL/functional_base.h>

namespace CGAL {

template <class Point_3>
class Less_xy_plane_yz_2 
{
public:
   typedef bool         result_type;
   typedef Arity_tag<2> Arity;

   bool 
   operator()(const Point_3& p, const Point_3& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.y(), p.z(), q.y(), q.z()) == SMALLER;
   }
};

template <class Point_3>
class Less_yx_plane_yz_2 
{
public:
   typedef bool         result_type;
   typedef Arity_tag<2> Arity;

   bool 
   operator()(const Point_3& p, const Point_3& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.z(), p.y(), q.z(), q.y()) == SMALLER;
   }
};

template <class Point_3>
class Left_turn_plane_yz_2 
{
public:
   typedef bool         result_type;
   typedef Arity_tag<3> Arity;

   bool 
   operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
   { 
    return orientationC2(p.y(), p.z(), q.y(), q.z(), r.y(), r.z()) == LEFTTURN;
   }
};

template <class Point_3>
class Left_of_line_plane_yz_2
{
public:
   Left_of_line_plane_yz_2(const Point_3& a, const Point_3& b):
      p_a(a), p_b(b) 
   { }

   bool 
   operator()(const Point_3& c) const
   {
      return orientationC2(p_a.y(), p_a.z(), p_b.y(), p_b.z(), c.y(), c.z()) ==
             LEFTTURN; 
   }
private:
   Point_3 p_a;
   Point_3 p_b;
};

template <class Point_3>
class Less_dist_to_line_plane_yz_2
{
public:
   typedef bool         result_type;
   typedef Arity_tag<4> Arity;

   bool
   operator()(const Point_3& p, const Point_3& q,
              const Point_3& r, const Point_3& s) const
   {
      Comparison_result
         res = cmp_signed_dist_to_lineC2(p.y(), p.z(), q.y(), q.z(),
                                         r.y(), r.z(), s.y(), s.z());
      if ( res == LARGER )
         return false;
      else if ( res == SMALLER )
         return true;
      else
         return compare_lexicographically_xyC2(r.y(), r.z(), s.y(), s.z())
             == SMALLER;
   }
};

template <class Point_3>
class Less_rotate_ccw_plane_yz_2
{
public:
   typedef bool         result_type;
   typedef Arity_tag<3> Arity;

   bool
   operator()(const Point_3& r, const Point_3& p, const Point_3& q) const
   {
      Orientation orient =
               orientationC2(r.y(), r.z(), p.y(), p.z(), q.y(), q.z());
      if ( orient ==  LEFTTURN )
         return true;
      else if ( orient == RIGHTTURN )
         return false;
      else
      {
         if (p.y() == r.y() && p.z() == r.z()) return false;
         if (q.y() == r.y() && q.z() == r.z()) return true;
         if (p.y() == q.y() && p.z() == q.z()) return false;
         return
            collinear_are_ordered_along_lineC2(r.y(), r.z(),
                                               q.y(), q.z(), p.y(), p.z());
      }
   }
};



template <class Point_3>
class Convex_hull_projective_yz_traits_2 
{
public:
    typedef Point_3                            Point_2;
    typedef Less_xy_plane_yz_2<Point_3>        Less_xy_2;
    typedef Less_yx_plane_yz_2<Point_3>        Less_yx_2;
    typedef Left_turn_plane_yz_2<Point_3>      Leftturn_2;
    typedef Less_dist_to_line_plane_yz_2<Point_3>
                                                Less_signed_distance_to_line_2;
    typedef Less_rotate_ccw_plane_yz_2<Point_3> Less_rotate_ccw_plane_2;

    Less_xy_2
    less_xy_2_object() const
    {  return Less_xy_2(); }

    Less_yx_2
    less_yx_2_object() const
    {  return Less_yx_2(); }

    Leftturn_2
    leftturn_2_object() const
    {  return Leftturn_2(); }

    Less_signed_distance_to_line_2
    less_signed_distance_to_line_2_object() const
    {  return Less_signed_distance_to_line_2(); }

    Less_rotate_ccw_plane_2
    less_rotate_ccw_plane_2_object() const
    {  return Less_rotate_ccw_plane_2(); }
};

}
#endif // CONVEX_HULL_PROJECTIVE_YZ_TRAITS_2_H
