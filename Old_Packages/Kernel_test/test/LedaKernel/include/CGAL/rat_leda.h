// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/rat_leda.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_RAT_LEDA_H
#define CGAL_RAT_LEDA_H

#include <CGAL/user_classes.h>
#include <CGAL/basic_classes.h>
#include <LEDA/rat_point.h>
#include <CGAL/Homogeneous.h>

#define CGAL_REP_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE


class use_rat_leda_kernel;
CGAL_TEMPLATE_NULL class Point_2<use_rat_leda_kernel>;

class use_rat_leda_kernel : public Homogeneous_base< use_rat_leda_kernel, leda_integer, leda_rational>
{
  public:
    typedef  use_rat_leda_kernel  R;

    typedef  leda_integer      RT;
    typedef  leda_rational     FT;

    typedef  leda_rat_point    Point_2_base;
    typedef  leda_rat_vector   Vector_2_base;

    typedef  CGAL::Point_2<R>  Point_2;
    typedef  CGAL::Vector_2<R> Vector_2;
};


CGAL_TEMPLATE_NULL
class Point_2<use_rat_leda_kernel> : public leda_rat_point
{
  public:

    typedef use_rat_leda_kernel      R;
    typedef use_rat_leda_kernel::RT  RT;
    typedef use_rat_leda_kernel::FT  FT;

    Point_2<use_rat_leda_kernel>()
    {}

    Point_2<use_rat_leda_kernel>(const Origin &o)
    {}

    Point_2<use_rat_leda_kernel>(const leda_rat_point& p)
      : leda_rat_point(p)
    {}

    Point_2<use_rat_leda_kernel>(const RT &hx, const RT &hy)
      : leda_rat_point(hx, hy)
    {}

    Point_2<use_rat_leda_kernel>(const RT &hx, const RT &hy, const RT &hw)
      : leda_rat_point(hx, hy, hw)
    {}
};

inline
Comparison_result
compare_distance_to_point( const leda_rat_point& p,
                   const leda_rat_point& q,
                   const leda_rat_point& r)
{
  if ( q == r) return EQUAL;
  return ( p.sqr_dist(q) > p.sqr_dist(r) ) ? LARGER : SMALLER;
}

inline
bool
has_larger_distance_to_point( const leda_rat_point& p,
                          const leda_rat_point& q,
                          const leda_rat_point& r)
{ return (compare_distance_to_point(p,q,r) == LARGER); }

inline
bool
has_smaller_distance_to_point( const leda_rat_point& p,
                           const leda_rat_point& q,
                           const leda_rat_point& r)
{ return (compare_distance_to_point(p,q,r) == SMALLER); }

inline
bool
has_smaller_signed_distance_to_line( const leda_rat_point& p,
                                     const leda_rat_point& q,
                                     const leda_rat_point& r,
                                     const leda_rat_point& s)
{ return (cmp_signed_dist(p,q,r,s) == SMALLER); }

inline
bool
has_larger_signed_distance_to_line( const leda_rat_point& p,
                                     const leda_rat_point& q,
                                     const leda_rat_point& r,
                                     const leda_rat_point& s)
{ return (cmp_signed_dist(p,q,r,s) == LARGER); }

inline
Comparison_result
compare_signed_distance_to_line( const leda_rat_point& p,
                                     const leda_rat_point& q,
                                     const leda_rat_point& r,
                                     const leda_rat_point& s)
{ return Comparison_result(cmp_signed_dist(p,q,r,s)); }

inline
bool 
left_turn(const leda_rat_point& p, const leda_rat_point& q, 
          const leda_rat_point& r)
{
   return ::left_turn(p,q,r);
}

inline
bool 
right_turn(const leda_rat_point& p, const leda_rat_point& q, 
          const leda_rat_point& r)
{
   return ::right_turn(p,q,r);
}

inline
bool
lexicographically_yx_larger_or_equal(const leda_rat_point& p,
                                     const leda_rat_point& q)
{
   if (p.ycoord() > q.ycoord()) return true;
   if (p.ycoord() < q.ycoord()) return false;
   return (p.xcoord() >= q.xcoord()); 
}

inline
bool
lexicographically_xy_larger_or_equal(const leda_rat_point& p,
                                     const leda_rat_point& q)
{
   if (p.xcoord() > q.xcoord()) return true;
   if (p.xcoord() < q.xcoord()) return false;
   return (p.ycoord() >= q.ycoord()); 
}

inline
Comparison_result
compare_yx(const leda_rat_point& p, const leda_rat_point& q)
{
   if (p.ycoord() > q.ycoord()) return CGAL::LARGER;
   if (p.ycoord() < q.ycoord()) return CGAL::SMALLER;
   if (p.xcoord() > q.xcoord()) return CGAL::LARGER; 
   if (p.xcoord() < q.xcoord()) return CGAL::SMALLER; 
   return CGAL::EQUAL;
}

inline
Comparison_result
compare_xy(const leda_rat_point& p, const leda_rat_point& q)
{
   if (p.xcoord() > q.xcoord()) return CGAL::LARGER; 
   if (p.xcoord() < q.xcoord()) return CGAL::SMALLER; 
   if (p.ycoord() > q.ycoord()) return CGAL::LARGER;
   if (p.ycoord() < q.ycoord()) return CGAL::SMALLER;
   return CGAL::EQUAL;
}
CGAL_END_NAMESPACE


#endif // CGAL_RAT_LEDA_H

