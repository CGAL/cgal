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
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>

#ifndef CGAL_IN_SMALLEST_ORTHOGONALCIRCLE_FTC2_H 
#define CGAL_IN_SMALLEST_ORTHOGONALCIRCLE_FTC2_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>

//-------------------------------------------------------------------
namespace CGAL {
//-------------------------------------------------------------------

template< class FT >
CGAL_MEDIUM_INLINE
Bounded_side
in_smallest_orthogonalcircleC2(const FT &px, const FT &py, const FT  &pw,
			       const FT &qx, const FT &qy, const FT  &qw,  
			       const FT &tx, const FT &ty, const FT  &tw)
{
  FT dpx = px-qx;
  FT dpy = py-qy;
  FT dtx = tx-qx;
  FT dty = ty-qy;
  FT dpz = CGAL_NTS square(dpx)+CGAL_NTS square(dpy);
 
  return Bounded_side 
    (CGAL_NTS sign(-(CGAL_NTS square(dtx)+CGAL_NTS square(dty)-tw+qw)*dpz
		   +(dpz-pw+qw)*(dpx*dtx+dpy*dty)));
}

//-------------------------------------------------------------------
} //namespace CGAL
//-------------------------------------------------------------------

#endif //CGAL_IN_SMALLEST_ORTHOGONALCIRCLEC2_H
