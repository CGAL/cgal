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

#ifndef CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALCIRCLE_FTC2_H 
#define CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALCIRCLE_FTC2_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>

//-------------------------------------------------------------------
namespace CGAL {
//-------------------------------------------------------------------

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_orthogonalcircleC2(
  const FT &px, const FT &py, const FT  &pw,
  const FT &qx, const FT &qy, const FT  &qw,  
  const FT &rx, const FT &ry, const FT  &rw)

{
  FT FT4(4);
  FT dpx = px-rx;
  FT dpy = py-ry;
  FT dqx = qx-rx;
  FT dqy = qy-ry;
  FT dpp = CGAL_NTS square(dpx)+CGAL_NTS square(dpy)-pw+rw;
  FT dqq = CGAL_NTS square(dqx)+CGAL_NTS square(dqy)-qw+rw;

  FT det0 = determinant(dpx, dpy, dqx, dqy);
  
  FT det1 = determinant(dpp, dpy, dqq, dqy);

  FT det2 = determinant(dpx, dpp, dqx, dqq);

  return 
    (CGAL_NTS square(det1)+CGAL_NTS square(det2))/
                                  (FT4*CGAL_NTS square(det0)) - rw;
}

template< class FT >
CGAL_MEDIUM_INLINE
FT
squared_radius_smallest_orthogonalcircleC2(
  const FT &px, const FT &py, const FT  &pw,
  const FT &qx, const FT &qy, const FT  &qw)

{
  FT FT4(4);
  FT dpz = CGAL_NTS square(px-qx)+CGAL_NTS square(py-qy);

  return (CGAL_NTS square(dpz-pw+qw)/(FT4*dpz)-qw);
}

//-------------------------------------------------------------------
} //namespace CGAL
//-------------------------------------------------------------------

#endif //CGAL_SQUARED_RADIUS_SMALLEST_ORTHOGONALCIRCLE_ftC2_H
