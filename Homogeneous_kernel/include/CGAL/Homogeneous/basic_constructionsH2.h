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
// SPDX-License-Identifier: LGPL-3.0+
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


template < class RT, class We>
void
weighted_circumcenterH2( const RT &phx, const RT &phy, const RT &phw, 
			 const We &pwt,
			 const RT &qhx, const RT &qhy, const RT &qhw, 
			 const We &qwt,
			 const RT &rhx, const RT &rhy, const RT &rhw, 
			 const We &rwt,
			 RT &vvx, RT &vvy, RT &vvw )
{

//   RT qx_px = ( qhx*qhw*phw*phw - phx*phw*qhw*qhw );
//   RT qy_py = ( qhy*qhw*phw*phw - phy*phw*qhw*qhw );
//   RT rx_px = ( rhx*rhw*phw*phw - phx*phw*rhw*rhw );
//   RT ry_py = ( rhy*rhw*phw*phw - phy*phw*rhw*rhw );

//    //intersection of the two radical axis of (qp) and (rp)
//   RT px2_py2_rx2_ry_2 =
//     phx*phx*rhw*rhw + phy*phy*rhw*rhw - rhx*rhx*phw*phw -
//     rhy*rhy*phw*phw - RT(pwt*pwt) + RT(rwt*rwt);
//   RT px2_py2_qx2_qy_2 =
//     phx*phx*qhw*qhw + phy*phy*qhw*qhw - qhx*qhx*phw*phw -
//     qhy*qhy*phw*phw - RT(pwt*pwt) + RT(qwt*qwt);
  
//   vvx = qy_py * px2_py2_rx2_ry_2 - ry_py * px2_py2_qx2_qy_2;
//   vvy = rx_px * px2_py2_qx2_qy_2 - qx_px * px2_py2_rx2_ry_2;
//   vvw = RT(2) * ( qx_px * ry_py - rx_px * qy_py );

  RT a1, b1, c1;
  RT a2, b2, c2;
  radical_axisH2(phx,phy,phw,pwt,qhx,qhy,qhw,qwt,a1,b1,c1);
  radical_axisH2(phx,phy,phw,pwt,rhx,rhy,rhw,rwt,a2,b2,c2);
  vvx = b1*c2 - c1*b2;
  vvy = c1*a2 - c2*a1;
  vvw = a1*b2 - a2*b1;
}




template < class RT , class FT>
void
radical_axisH2(const RT &phx, const RT &phy, const RT &phw, const FT &pwt,
	       const RT &qhx, const RT &qhy, const RT &qhw, const FT &qwt,
	       RT &a, RT &b, RT &c )
{
  Rational_traits<FT> rt;
  RT npwt = rt.numerator(pwt);
  RT dpwt = rt.denominator(pwt);
  RT nqwt = rt.numerator(qwt);
  RT dqwt = rt.denominator(qwt);

  a =  RT(2) * ( phx*phw*qhw*qhw - qhx*qhw*phw*phw );
  b =  RT(2) * ( phy*phw*qhw*qhw - qhy*qhw*phw*phw );
  c = (- phx*phx*qhw*qhw - phy*phy*qhw*qhw 
       + qhx*qhx*phw*phw + qhy*qhy*phw*phw) * dpwt * dqwt 
      + npwt*dqwt*phw*phw*qhw*qhw - nqwt*dpwt*phw*phw*qhw*qhw; 

}



} //namespace CGAL

#endif // CGAL_BASIC_CONSTRUCTIONSH2_H
