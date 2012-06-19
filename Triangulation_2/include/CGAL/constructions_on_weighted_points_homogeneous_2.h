// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Julia Flototto, Mariette Yvinec


#ifndef CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_HOMOGENEOUS_2_H
#define CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_HOMOGENEOUS_2_H

namespace CGAL {

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


template < class RT , class We>
void
radical_axisH2(const RT &phx, const RT &phy, const RT &phw, const We &pwt,
	       const RT &qhx, const RT &qhy, const RT &qhw, const We &qwt,
	       RT &a, RT &b, RT &c )
{
//   a = RT(2) * ( qhx*qhw*phw*phw - phx*phw*qhw*qhw );
//   b = RT(2) * ( qhy*qhw*phw*phw - phy*phw*qhw*qhw );
//   c = phx*phx*qhw*qhw + phy*phy*qhw*qhw - qhx*qhx*phw*phw 
//     - qhy*qhy*phw*phw- RT(pwt*pwt) + RT(qwt*qwt); 

  a =  RT(2) * ( phx*phw*qhw*qhw - qhx*qhw*phw*phw );
  b =  RT(2) * ( phy*phw*qhw*qhw - qhy*qhw*phw*phw );
  c = - phx*phx*qhw*qhw - phy*phy*qhw*qhw 
      + qhx*qhx*phw*phw + qhy*qhy*phw*phw 
      + RT(pwt)*phw*phw*qhw*qhw - RT(qwt)*phw*phw*qhw*qhw; 

}




} //namespace CGAL
#endif //CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_HOMOGENEOUS_2_H
