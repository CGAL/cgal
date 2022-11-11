// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra


#ifndef CGAL_BASIC_CONSTRUCTIONSH3_H
#define CGAL_BASIC_CONSTRUCTIONSH3_H

#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/Cartesian/Plane_3.h>


namespace CGAL {

template <class R>
typename R::Point_3
_projection(const typename R::Point_3& p, const PlaneH3<R>& pl)
{
  typedef typename R::RT RT;
  if ( pl.has_on(p) ) return p;

  RT A = pl.a();
  RT B = pl.b();
  RT C = pl.c();
  RT D = pl.d();
  RT phx = p.hx();
  RT phy = p.hy();
  RT phz = p.hz();
  RT phw = p.hw();

  RT num = A * phx  +  B * phy  +  C * phz  +  D * phw;
  RT den = A * A    +  B * B    +  C * C;

  return typename R::Point_3( num * A  -  den * phx,
                              num * B  -  den * phy,
                              num * C  -  den * phz,
                             -den );
}

template <class R>
typename R::Point_3
gp_linear_intersection(const PlaneH3<R> &f,
                       const PlaneH3<R> &g,
                       const PlaneH3<R> &h)
{
  typedef typename R::RT RT;
  return typename R::Point_3(
                  determinant<RT>(-f.d(), f.b(), f.c(),
                                        -g.d(), g.b(), g.c(),
                                        -h.d(), h.b(), h.c()),
                  determinant<RT>( f.a(),-f.d(), f.c(),
                                         g.a(),-g.d(), g.c(),
                                         h.a(),-h.d(), h.c()),
                  determinant<RT>( f.a(), f.b(),-f.d(),
                                         g.a(), g.b(),-g.d(),
                                         h.a(), h.b(),-h.d()),
                  determinant<RT>( f.a(), f.b(), f.c(),
                                         g.a(), g.b(), g.c(),
                                         h.a(), h.b(), h.c()));
}

template <class R>
CGAL_KERNEL_INLINE
typename R::FT
squared_distance( PointH3<R> const& p, PointH3<R> const& q)
{ return (p-q)*(p-q); }

} //namespace CGAL

#endif // CGAL_BASIC_CONSTRUCTIONSH3_H
