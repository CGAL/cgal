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
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_ORIENTATION_PREDICATESH3_H
#define CGAL_ORIENTATION_PREDICATESH3_H

#include <CGAL/Homogeneous/PointH3.h>
#include <CGAL/Homogeneous/predicates_on_rtH2.h>
#include <CGAL/predicates/sign_of_determinant.h>

CGAL_BEGIN_NAMESPACE

template < class R>
CGAL_KERNEL_INLINE
Orientation
orientation( const PointH3<R>& p,
             const PointH3<R>& q,
             const PointH3<R>& r,
             const PointH3<R>& s)
{
  // Two rows are switched, because of the homogeneous column.
  return (Orientation) sign_of_determinant4x4( p.hx(), p.hy(), p.hz(), p.hw(),
                                               r.hx(), r.hy(), r.hz(), r.hw(),
                                               q.hx(), q.hy(), q.hz(), q.hw(),
                                               s.hx(), s.hy(), s.hz(), s.hw());
}

template < class R>
inline
bool
are_positive_oriented( const PointH3<R>& p,
                       const PointH3<R>& q,
                       const PointH3<R>& r,
                       const PointH3<R>& s)
{ return (orientation(p,q,r,s) == POSITIVE); }

template < class R>
inline
bool
are_negative_oriented( const PointH3<R>& p,
                       const PointH3<R>& q,
                       const PointH3<R>& r,
                       const PointH3<R>& s)
{ return (orientation(p,q,r,s) == NEGATIVE); }

template < class R>
inline
bool
coplanar( const PointH3<R>& p,
          const PointH3<R>& q,
          const PointH3<R>& r,
          const PointH3<R>& s)
{ return (orientation(p,q,r,s) == COPLANAR); }


template <class R>
Orientation
coplanar_orientation(const PointH3<R>& p,
                     const PointH3<R>& q,
                     const PointH3<R>& r,
                     const PointH3<R>& s)
// p,q,r,s supposed to be coplanar
// p, q, r supposed to be non collinear
// tests whether s is on the same side of p,q as r
// returns :
// COLLINEAR if pqs collinear
// POSITIVE if pqr and pqs have the same orientation
// NEGATIVE if pqr and pqs have opposite orientations
{
  CGAL_kernel_precondition( coplanar( p, q, r, s));
  // p, q, r define a plane P:
  CGAL_kernel_precondition( !collinear( p, q, r));
  // compute orientation of p,q,s in this plane P:
  Orientation save;
  if ( (save = orientationH2( p.hy(), p.hz(), p.hw(),
                              q.hy(), q.hz(), q.hw(),
                              r.hy(), r.hz(), r.hw())) != COLLINEAR)
  { return
      static_cast<Orientation>(
        static_cast<int>( save)
      * static_cast<int>( orientationH2( p.hy(), p.hz(), p.hw(),
                                         q.hy(), q.hz(), q.hw(),
                                         s.hy(), s.hz(), s.hw())) );
  }
  if ( (save = orientationH2( p.hx(), p.hz(), p.hw(),
                              q.hx(), q.hz(), q.hw(),
                              r.hx(), r.hz(), r.hw())) != COLLINEAR)
  { return
      static_cast<Orientation>(
        static_cast<int>( save)
      * static_cast<int>( orientationH2( p.hx(), p.hz(), p.hw(),
                                         q.hx(), q.hz(), q.hw(),
                                         s.hx(), s.hz(), s.hw())) );
  }
  if ( (save = orientationH2( p.hx(), p.hy(), p.hw(),
                              q.hx(), q.hy(), q.hw(),
                              r.hx(), r.hy(), r.hw())) != COLLINEAR)
  { return
      static_cast<Orientation>(
        static_cast<int>( save)
      * static_cast<int>( orientationH2( p.hx(), p.hy(), p.hw(),
                                         q.hx(), q.hy(), q.hw(),
                                         s.hx(), s.hy(), s.hw())) );
  }
  CGAL_kernel_assertion( false);
  return COLLINEAR;
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
Orientation
coplanar_orientation(const PointH3<R>& p,
                     const PointH3<R>& q,
                     const PointH3<R>& r)
// Returns an Orientation which is coherent for all (p,q,r) chosen in a same
// plane.
{
  Orientation oxy_pqr = orientationH2(p.hx(), p.hy(), p.hw(),
	                              q.hx(), q.hy(), q.hw(),
				      r.hx(), r.hy(), r.hw());
  if (oxy_pqr != COLLINEAR)
      return oxy_pqr;

  Orientation oyz_pqr = orientationH2(p.hy(), p.hz(), p.hw(),
	                              q.hy(), q.hz(), q.hw(),
				      r.hy(), r.hz(), r.hw());
  if (oyz_pqr != COLLINEAR)
      return oyz_pqr;

  return orientationH2(p.hx(), p.hz(), p.hw(),
	               q.hx(), q.hz(), q.hw(),
		       r.hx(), r.hz(), r.hw());
}

CGAL_END_NAMESPACE

#endif // CGAL_ORIENTATION_PREDICATESH3_H
