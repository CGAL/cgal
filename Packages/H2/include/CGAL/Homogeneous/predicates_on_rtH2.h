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
 

#ifndef CGAL_PREDICATES_ON_RTH2_H
#define CGAL_PREDICATES_ON_RTH2_H

CGAL_BEGIN_NAMESPACE

template <class RT>
CGAL_KERNEL_INLINE
Orientation
orientationH2( const RT& phx, const RT& phy, const RT& phw,
               const RT& qhx, const RT& qhy, const RT& qhw,
               const RT& rhx, const RT& rhy, const RT& rhw )
{
  const RT  RT0 = RT(0);
  
  // | A B |
  // | C D |
  
  RT  A = phx*rhw - phw*rhx;
  RT  B = phy*rhw - phw*rhy;
  RT  C = qhx*rhw - qhw*rhx;
  RT  D = qhy*rhw - qhw*rhy;
  
  RT  det =  A*D - B*C;
  
  
  if (det < RT0  )
  {
      return CLOCKWISE;
  }
  else
  {
      return (RT0 < det) ? COUNTERCLOCKWISE : COLLINEAR;
  }
}

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_ON_RTH2_H
