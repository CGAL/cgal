// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra
 
#ifndef CGAL_PREDICATES_ON_LINESH2_H
#define CGAL_PREDICATES_ON_LINESH2_H

#include <CGAL/Homogeneous/PointH2.h>
#include <CGAL/Homogeneous/LineH2.h>
#include <CGAL/Homogeneous/predicates_on_pointsH2.h>
#include <CGAL/Homogeneous/basic_constructionsH2.h>

CGAL_BEGIN_NAMESPACE

template <class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_x(const PointH2<R>& p,
          const LineH2<R>& l1,
          const LineH2<R>& l2)
{
  PointH2<R> ip = gp_linear_intersection( l1, l2 );
  return compare_x( p, ip );
}

template <class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_x(const LineH2<R>& l1,
          const LineH2<R>& l2,
          const LineH2<R>& h1,
          const LineH2<R>& h2)
{
  PointH2<R> lip = gp_linear_intersection( l1, l2 );
  PointH2<R> hip = gp_linear_intersection( h1, h2 );
  return compare_x( lip, hip );
}

template < class R >
inline
Comparison_result
compare_x(const LineH2<R> &l,
	  const LineH2<R> &h1,
	  const LineH2<R> &h2)
{
    return compare_x(l, h1, l, h2);
}

template <class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_y(const PointH2<R>& p,
          const LineH2<R>& l1,
          const LineH2<R>& l2)
{
  PointH2<R> ip = gp_linear_intersection( l1, l2 );
  return compare_y( p, ip );
}

template <class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_y(const LineH2<R>& l1,
          const LineH2<R>& l2,
          const LineH2<R>& h1,
          const LineH2<R>& h2)
{
  PointH2<R> lip = gp_linear_intersection( l1, l2 );
  PointH2<R> hip = gp_linear_intersection( h1, h2 );
  return compare_y( lip, hip );
}

template < class R >
inline
Comparison_result
compare_y(const LineH2<R> &l,
	  const LineH2<R> &h1,
	  const LineH2<R> &h2)
{
    return compare_y(l, h1, l, h2);
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_y_at_x(const PointH2<R>& p,
               const LineH2<R>& h)
{
  typedef typename R::RT RT;
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

template <class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_y_at_x(const PointH2<R>& p,
               const LineH2<R>& h1,
               const LineH2<R>& h2)
{ return CGAL_NTS compare(h1.y_at_x( p.x() ), h2.y_at_x( p.x() )); }

template <class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_y_at_x(const LineH2<R>& l1,
               const LineH2<R>& l2,
               const LineH2<R>& h)
{ return compare_y_at_x( gp_linear_intersection( l1, l2 ), h); }

template <class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_y_at_x(const LineH2<R>& l1,
               const LineH2<R>& l2,
               const LineH2<R>& h1,
               const LineH2<R>& h2)
{ return compare_y_at_x( gp_linear_intersection( l1, l2 ), h1, h2 ); }

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_x_at_y(const PointH2<R>& p,
               const LineH2<R>& h)
{
  typedef typename R::RT RT;
  CGAL_kernel_precondition( ! h.is_horizontal() );
  Oriented_side ors = h.oriented_side( p );
  if ( h.a() < RT(0) )
  {
      ors = opposite( ors );
  }
  if ( ors == ON_POSITIVE_SIDE )
  {
      return LARGER;
  }
  return ( ors == ON_NEGATIVE_SIDE ) ? SMALLER : EQUAL;
}

template <class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_x_at_y(const PointH2<R>& p,
               const LineH2<R>& h1,
               const LineH2<R>& h2)
{ return CGAL_NTS compare(h1.x_at_y( p.y() ), h2.x_at_y( p.y() )); }


template <class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_x_at_y(const LineH2<R>& l1,
               const LineH2<R>& l2,
               const LineH2<R>& h)
{ return compare_x_at_y( gp_linear_intersection( l1, l2 ), h); }

template <class R>
CGAL_KERNEL_INLINE
Comparison_result
compare_x_at_y(const LineH2<R>& l1,
               const LineH2<R>& l2,
               const LineH2<R>& h1,
               const LineH2<R>& h2)
{ return compare_x_at_y( gp_linear_intersection( l1, l2 ), h1, h2 ); }

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_slopes(const LineH2<R>& l1, const LineH2<R>& l2)
{
   typedef typename R::RT RT;
   if (l1.is_horizontal())
     return l2.is_vertical() ? 
         SMALLER : Comparison_result(CGAL_NTS sign<RT>(l2.a() * l2.b()));
   if (l2.is_horizontal()) 
     return l1.is_vertical() ? 
         LARGER : Comparison_result(- CGAL_NTS sign<RT>(l1.a() * l1.b()));
   if (l1.is_vertical()) return l2.is_vertical() ? EQUAL : LARGER;
   if (l2.is_vertical()) return SMALLER;
   int l1_sign = CGAL_NTS sign<RT>(-l1.a() * l1.b());
   int l2_sign = CGAL_NTS sign<RT>(-l2.a() * l2.b());

   if (l1_sign < l2_sign) return SMALLER;
   if (l1_sign > l2_sign) return LARGER;

   if (l1_sign > 0)
     return CGAL_NTS compare( CGAL_NTS abs<RT>(l1.a() * l2.b()),
                              CGAL_NTS abs<RT>(l2.a() * l1.b()) );

   return CGAL_NTS compare( CGAL_NTS abs<RT>(l2.a() * l1.b()),
                            CGAL_NTS abs<RT>(l1.a() * l2.b()) );
}

CGAL_END_NAMESPACE

#endif  // CGAL_PREDICATES_ON_LINESH2_H
