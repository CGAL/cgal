// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
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

CGAL_END_NAMESPACE

#endif // CGAL_ORIENTATION_PREDICATESH3_H
