// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 21
//
// file          : 
// package       : Convex_hull_3 (2.6)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : chull_traits.lw
// revision      : 2.3  
// revision_date : 01 Feb 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_CHULL_TRAITS_2_H
#define CGAL_CHULL_TRAITS_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <LEDA/array.h>
#include <LEDA/point.h>

CGAL_BEGIN_NAMESPACE
template <class _R>
class chull_traits_2
{
  public:
    typedef _R                                R;
    typedef typename _R::RT                   RT;
    typedef Point_2<_R>                       POINT;
    typedef POINT                             IPOINT;
    typedef Line_2<_R>                        PLANE;
    typedef Vector_2<_R>                      VECTOR;
    typedef Ray_2<_R>                         IRAY;

  static int
    side(const Line_2<_R>& pl, 
         const Point_2<_R>& p)
    { return (int)pl.oriented_side(p); }

  static bool
    affinely_independent(const leda_array<Point_2<_R> >& A)
    {
      int a = A.size();
      if (a > 3)
          return false;
      if (a == 3)
          return !collinear( A[0], A[1], A[2] );
      if (a == 2)
          return (A[0] != A[1] );  
      return true;
    }
    
  static bool
    contained_in_simplex(const leda_array<Point_2<_R> >& A,
                         const Point_2<_R> p)
    {
      int a = A.size();
      CGAL_assertion( a <= 3 );
      if (a == 3)
      {
          Triangle_2<_R> t( A[0], A[1], A[2] );
          return !t.has_on_unbounded_side(p);
      }
      if (a == 2)
      {
          Segment_2<_R> s( A[0], A[1] );
          return s.has_on(p);
      }
      if (a == 1)
      {
          return ( A[0] == p);
      }
    }
    
  static bool
    contained_in_affine_hull(const leda_array<Point_2<_R> >& A,
                             const Point_2<_R> p)
    {
      int a = A.size();
      CGAL_assertion( affinely_independent(A) );
      if (a >  2)
          return true;
      if (a == 2)
          return collinear( p, A[0], A[1] );
      if (a == 1)
          return ( p == A[0] );
      return false;
    }
    
  static Vector_2<_R>
    normal(const Line_2<_R>& l)
    { return Vector_2<_R>( l.a(), l.b(), l.c() ); }

  static Vector_2<_R>
    to_vector(const Point_2<_R>& p)
    { return p - ORIGIN; }

  static Vector_2<_R>
    zero_vector(int i)
    { 
      CGAL_assertion( i == 2 );
      return Vector_2<_R>(RT(0), RT(0) ); 
    }

  static Point_2<_R>
    to_ipoint(const Vector_2<_R>& v)
    { return ORIGIN + v; }

  static Line_2<_R>
    hyperplane_through_points(const leda_array<Point_2<_R> >& A,
                              const Point_2<_R>& p,
                              int side)
    {
      int a = A.size();
      CGAL_assertion( a <= 2 );
      Line_2<_R> l;
      if (a == 2)
      {
          l = Line_2<_R>( A[0], A[1] );
      }
      if (a == 1)
      {
          l = Line_2<_R>( A[0], A[0] + (ORIGIN - p) );
      }
      if (side != 0)
      {
          if ( (int)l.oriented_side(p) !=  side ) { l = l.opposite(); }
      }
      return l;
    }

  static bool
    intersection(const Line_2<_R>& l, 
                 const Ray_2<_R>& r,
                 Point_2<_R>& ip)
    {
      Object res = intersection(l, r);
      return assign(ip, res);
    }
  static int
    orientation( const Point_2<_R>& p1, 
                 const Point_2<_R>& p2,
                 const Point_2<_R>& p3 )
    { return (int)CGAL::orientation(p1,p2,p3); }

  static leda_point
    to_d2_point( const Point_2<_R>& p )
    {
      return leda_point( CGAL::to_double( p.x() ),
                         CGAL::to_double( p.y() ) );
    }
};
CGAL_END_NAMESPACE

#endif // CGAL_CHULL_TRAITS_2_H
