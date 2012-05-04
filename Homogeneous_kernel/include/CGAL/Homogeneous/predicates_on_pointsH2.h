// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_PREDICATES_ON_POINTSH2_H
#define CGAL_PREDICATES_ON_POINTSH2_H

#include <CGAL/Homogeneous/PointH2.h>

namespace CGAL {

template < class R>
CGAL_KERNEL_INLINE
bool
equal_xy(const PointH2<R>& p,
         const PointH2<R>& q)
{
  typedef typename R::RT RT;

  // Using these references allows to spare calls to [pq].hw().
  const RT& phw = p.hw();
  const RT& qhw = q.hw();

  return (p.hx()*qhw == q.hx()*phw) && (p.hy()*qhw == q.hy()*phw);
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
_where_wrt_L_wedge( const PointH2<R>& p, const PointH2<R>& q )
{
  Sign xs = CGAL_NTS sign( q.hx()*p.hw() - p.hx()*q.hw() );  // sign( qx - px )
  Sign ys = CGAL_NTS sign( q.hy()*p.hw() - p.hy()*q.hw() );  // sign( qy - py )

  if (( xs == NEGATIVE ) || ( ys == NEGATIVE ))
      return ON_NEGATIVE_SIDE;
  if (( xs == POSITIVE ) && ( ys == POSITIVE ))
      return ON_POSITIVE_SIDE;
  return ON_ORIENTED_BOUNDARY;
}

#if 0
// Unused, undocumented, un-functorized.
template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Comparison_result
compare_deltax_deltay(const PointH2<R>& p,
                      const PointH2<R>& q,
                      const PointH2<R>& r,
                      const PointH2<R>& s)
{
  return CGAL_NTS compare(
                  CGAL_NTS abs(p.hx()*q.hw() - q.hx()*p.hw()) * r.hw()*s.hw(),
                  CGAL_NTS abs(r.hy()*s.hw() - s.hy()*r.hw()) * p.hw()*q.hw());
}
#endif

} //namespace CGAL

#endif // CGAL_PREDICATES_ON_POINTSH2_H
