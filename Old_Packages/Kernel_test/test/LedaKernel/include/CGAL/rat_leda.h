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
// source        : test_kernel_programs.fw
// file          : include/CGAL/rat_leda.h
// revision      : 3.8
// revision_date : 08 Oct 2000 
// author(s)     : Stefan Schirra
//
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de> 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_RAT_LEDA_H
#define CGAL_RAT_LEDA_H

#include <CGAL/user_classes.h>
#include <CGAL/basic_classes.h>
#include <LEDA/rat_point.h>

#define CGAL_REP_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE


class use_rat_leda_kernel;
CGAL_TEMPLATE_NULL class Point_2<use_rat_leda_kernel>;

class use_rat_leda_kernel
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
cmp_dist_to_point( const leda_rat_point& p,
                   const leda_rat_point& q,
                   const leda_rat_point& r)
{
  if ( q == r) return EQUAL;
  return ( p.sqr_dist(q) > p.sqr_dist(r) ) ? LARGER : SMALLER;
}

inline
bool
has_larger_dist_to_point( const leda_rat_point& p,
                          const leda_rat_point& q,
                          const leda_rat_point& r)
{ return (cmp_dist_to_point(p,q,r) == LARGER); }

inline
bool
has_smaller_dist_to_point( const leda_rat_point& p,
                           const leda_rat_point& q,
                           const leda_rat_point& r)
{ return (cmp_dist_to_point(p,q,r) == SMALLER); }


CGAL_END_NAMESPACE


#endif // CGAL_RAT_LEDA_H

