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
// Author(s)     : Sven Schoenherr
//                 Stefan Schirra
 

#ifndef CGAL_BASIC_CONSTRUCTIONSH2_H
#define CGAL_BASIC_CONSTRUCTIONSH2_H

#include <CGAL/Homogeneous/PointH2.h>
#include <CGAL/Homogeneous/LineH2.h>

namespace CGAL {

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
typename R::Point_2
gp_linear_intersection(const LineH2<R>& l1, const LineH2<R>& l2)
{
  return typename R::Point_2( l1.b()*l2.c() - l2.b()*l1.c(),
                              l2.a()*l1.c() - l1.a()*l2.c(),
                              l1.a()*l2.b() - l2.a()*l1.b() );
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
typename R::FT
squared_distance( const PointH2<R>& p, const PointH2<R>& q )
{
  typedef typename R::RT RT;
  typedef typename R::FT FT;

  const RT & phx = p.hx();
  const RT & phy = p.hy();
  const RT & phw = p.hw();
  const RT & qhx = q.hx();
  const RT & qhy = q.hy();
  const RT & qhw = q.hw();

  RT sq_dist_numerator =
          phx * phx * qhw * qhw
    - RT(2) * phx * qhx * phw * qhw
    +     qhx * qhx * phw * phw

    +     phy * phy * qhw * qhw
    - RT(2) * phy * qhy * phw * qhw
    +     qhy * qhy * phw * phw ;

  RT sq_dist_denominator = qhw * qhw * phw * phw ;

  return FT( sq_dist_numerator ) / FT( sq_dist_denominator );
}

} //namespace CGAL

#endif // CGAL_BASIC_CONSTRUCTIONSH2_H
