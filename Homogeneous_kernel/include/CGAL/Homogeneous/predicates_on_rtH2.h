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
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra
 

#ifndef CGAL_PREDICATES_ON_RTH2_H
#define CGAL_PREDICATES_ON_RTH2_H

#include <CGAL/predicates/sign_of_determinant.h>

CGAL_BEGIN_NAMESPACE

template <class RT>
inline
Orientation
orientationH2( const RT& phx, const RT& phy, const RT& phw,
               const RT& qhx, const RT& qhy, const RT& qhw,
               const RT& rhx, const RT& rhy, const RT& rhw )
{
  return sign_of_determinant3x3(phx,phy,phw,
                                qhx,qhy,qhw,
                                rhx,rhy,rhw);
}

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_ON_RTH2_H
