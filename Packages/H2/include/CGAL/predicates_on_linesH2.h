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
// release_date  : 2000, October 11
// 
// source        : predicates_on_linesH2.fw
// file          : predicates_on_linesH2.h
// package       : H2 (2.13)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.13
// revision_date : 11 Oct 2000 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_PREDICATES_ON_LINESH2_H
#define CGAL_PREDICATES_ON_LINESH2_H

#ifndef CGAL_POINTH2_H
#include <CGAL/PointH2.h>
#endif // CGAL_POINTH2_H
#ifndef CGAL_LINEH2_H
#include <CGAL/LineH2.h>
#endif // CGAL_LINEH2_H
#ifndef CGAL_PREDICATES_ON_POINTSH2_H
#include <CGAL/predicates_on_pointsH2.h>
#endif // CGAL_PREDICATES_ON_POINTSH2_H
#ifndef CGAL_BASIC_CONSTRUCTIONSH2_H
#include <CGAL/basic_constructionsH2.h>
#endif // CGAL_BASIC_CONSTRUCTIONSH2_H

CGAL_BEGIN_NAMESPACE


template <class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
compare_x(const PointH2<FT,RT>& p,
          const LineH2<FT,RT>& l1,
          const LineH2<FT,RT>& l2)
{
  PointH2<FT,RT> ip = gp_linear_intersection( l1, l2 );
  return compare_x( p, ip );
}


template <class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
compare_x(const LineH2<FT,RT>& l1,
          const LineH2<FT,RT>& l2,
          const LineH2<FT,RT>& h1,
          const LineH2<FT,RT>& h2)
{
  PointH2<FT,RT> lip = gp_linear_intersection( l1, l2 );
  PointH2<FT,RT> hip = gp_linear_intersection( h1, h2 );
  return compare_x( lip, hip );
}
template <class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
compare_y(const PointH2<FT,RT>& p,
          const LineH2<FT,RT>& l1,
          const LineH2<FT,RT>& l2)
{
  PointH2<FT,RT> ip = gp_linear_intersection( l1, l2 );
  return compare_y( p, ip );
}


template <class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
compare_y(const LineH2<FT,RT>& l1,
          const LineH2<FT,RT>& l2,
          const LineH2<FT,RT>& h1,
          const LineH2<FT,RT>& h2)
{
  PointH2<FT,RT> lip = gp_linear_intersection( l1, l2 );
  PointH2<FT,RT> hip = gp_linear_intersection( h1, h2 );
  return compare_y( lip, hip );
}
template <class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_y_at_x(const PointH2<FT,RT>& p,
               const LineH2<FT,RT>& h)
{
  CGAL_kernel_precondition( ! h.is_vertical() );
  Oriented_side ors = h.oriented_side( p );
  if ( h.b() < RT(0) )
  {
      ors = opposite( ors );
  }
  if ( ors == ON_POSITIVE_SIDE )
  {
      return LARGER;
  }
  return ( ors == ON_NEGATIVE_SIDE ) ? SMALLER : EQUAL;
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
compare_y_at_x(const PointH2<FT,RT>& p,
               const LineH2<FT,RT>& h1,
               const LineH2<FT,RT>& h2)
{ return CGAL_NTS compare(h1.y_at_x( p.x() ), h2.y_at_x( p.x() )); }


template <class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
compare_y_at_x(const LineH2<FT,RT>& l1,
               const LineH2<FT,RT>& l2,
               const LineH2<FT,RT>& h)
{ return compare_y_at_x( gp_linear_intersection( l1, l2 ), h); }


template <class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
compare_y_at_x(const LineH2<FT,RT>& l1,
               const LineH2<FT,RT>& l2,
               const LineH2<FT,RT>& h1,
               const LineH2<FT,RT>& h2)
{ return compare_y_at_x( gp_linear_intersection( l1, l2 ), h1, h2 ); }

template <class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_x_at_y(const PointH2<FT,RT>& p,
               const LineH2<FT,RT>& h)
{
  CGAL_kernel_precondition( ! h.is_horizontal() );
  Oriented_side ors = h.oriented_side( p );
  if ( h.b() < RT(0) )
  {
      ors = opposite( ors );
  }
  if ( ors == ON_POSITIVE_SIDE )
  {
      return LARGER;
  }
  return ( ors == ON_NEGATIVE_SIDE ) ? SMALLER : EQUAL;
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
compare_x_at_y(const PointH2<FT,RT>& p,
               const LineH2<FT,RT>& h1,
               const LineH2<FT,RT>& h2)
{ return CGAL_NTS compare(h1.x_at_y( p.y() ), h2.x_at_y( p.y() )); }


template <class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
compare_x_at_y(const LineH2<FT,RT>& l1,
               const LineH2<FT,RT>& l2,
               const LineH2<FT,RT>& h)
{ return compare_x_at_y( gp_linear_intersection( l1, l2 ), h); }


template <class FT, class RT>
CGAL_KERNEL_INLINE
Comparison_result
compare_x_at_y(const LineH2<FT,RT>& l1,
               const LineH2<FT,RT>& l2,
               const LineH2<FT,RT>& h1,
               const LineH2<FT,RT>& h2)
{ return compare_x_at_y( gp_linear_intersection( l1, l2 ), h1, h2 ); }

CGAL_END_NAMESPACE


#endif  // CGAL_PREDICATES_ON_LINESH2_H
