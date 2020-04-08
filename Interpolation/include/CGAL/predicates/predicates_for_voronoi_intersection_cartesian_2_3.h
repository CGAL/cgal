// Copyright (c) 2003 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Julia Floetotto

#ifndef CGAL_PREDICATES_FOR_VORONOI_INTERSECTION_CARTESIAN_2_3_H
#define CGAL_PREDICATES_FOR_VORONOI_INTERSECTION_CARTESIAN_2_3_H

#include <CGAL/license/Interpolation.h>


#include <CGAL/number_utils.h>
#include <CGAL/predicates/sign_of_determinant.h>

namespace CGAL {

template < class RT>
Oriented_side
side_of_plane_centered_sphere_translateC3(
    const RT &ax, const RT &ay, const RT &az,
    const RT &nx, const RT &ny, const RT &nz,
    const RT &qx, const RT &qy, const RT &qz,
    const RT &rx, const RT &ry, const RT &rz,
    const RT &tx, const RT &ty, const RT &tz)
{
  RT q2 = CGAL_NTS square(qx) + CGAL_NTS square(qy) + CGAL_NTS square(qz);
  RT r2 = CGAL_NTS square(rx) + CGAL_NTS square(ry) + CGAL_NTS square(rz);
  RT t2 = CGAL_NTS square(tx) + CGAL_NTS square(ty) + CGAL_NTS square(tz);
  RT na =nx*ax + ny*ay + nz*az;
  na *= RT(2.0);

  Sign num = sign_of_determinant(rx, ry, rz, r2,
                                 qx, qy, qz, q2,
                                 nx, ny, nz, na,
                                 tx, ty, tz, t2);

  //denumerator:
  Sign  den = sign_of_determinant(nx,ny,nz,
                                  qx,qy,qz,
                                  rx,ry,rz);
  CGAL_assertion(den != ZERO);

  return den * num;
}

template < class RT>
Oriented_side
side_of_plane_centered_sphereC3(const RT &ax, const RT &ay, const RT &az,
                                const RT &nx, const RT &ny, const RT &nz,
                                const RT &px, const RT &py, const RT &pz,
                                const RT &qx, const RT &qy, const RT &qz,
                                const RT &rx, const RT &ry, const RT &rz,
                                const RT &tx, const RT &ty, const RT &tz)
{
  // resolution of the system (where c denotes the sphere's center)
  //
  // ((a-c)*n)= 0                //the center lies in the plane (a, n)
  // (c-p)(c-p)=(c-q)(c-q)=(c-r)(c-r)   // p,q,r on a sphere(center c)
  //
  //  return:   sign( (c-p)(c-p) - (c-t)(c-t))
  //
  //method:
  // - tranlation of p to the origin.
  // - separate computation of det and norm of the expression

  return side_of_plane_centered_sphere_translateC3<RT>(ax-px, ay-py, az-pz,
                                                       nx, ny, nz,
                                                       qx-px, qy-py,qz-pz,
                                                       rx-px, ry-py,rz-pz,
                                                       tx-px, ty-py,tz-pz);
}

template < class RT>
Oriented_side
side_of_plane_centered_sphere_translateC3(
    const RT &ax, const RT &ay, const RT &az,
    const RT &nx, const RT &ny, const RT &nz,
    const RT &qx, const RT &qy, const RT &qz,
    const RT &rx, const RT &ry, const RT &rz)
{
  //first choice of n_ortho: (ny+nz, -nx, -nx)
  // if it is
  RT q2 = CGAL_NTS square(qx) + CGAL_NTS square(qy) + CGAL_NTS square(qz);
  RT r2 = CGAL_NTS square(rx) + CGAL_NTS square(ry) + CGAL_NTS square(rz);
  RT na =nx*ax + ny*ay + nz*az;
  na *= RT(2.0);

  Sign num = sign_of_determinant<RT>(qx, qy, qz, q2,
                                     ny, -nx, 0, 0,
                                     nx, ny, nz, na,
                                     rx, ry, rz, r2);
  //denumerator:
  Sign  den = sign_of_determinant<RT>(nx,ny,nz,
                                      ny,-nx, 0,
                                      qx,qy,qz);
  if (den==ZERO) {
    // bad choice: (ny,-nx,0) is coplanar with n,q.
    // by precondition: q and n may not be collinear
    // => the cross product q*n is orthogonal to q, n and not coplanar
    num = sign_of_determinant<RT>(qx, qy, qz, q2,
                                  ny*qz-nz*qy, nz*qx-nx*qz,nx*qy-ny*qx, 0,
                                  nx, ny, nz, na,
                                  rx, ry, rz, r2);
    den = sign_of_determinant<RT>(nx,ny,nz,
                                  ny*qz-nz*qy, nz*qx - nx*qz,nx*qy-ny*qx,
                                  qx,qy,qz);
  }
  CGAL_assertion(den != ZERO);
  return den * num;
}

template < class RT>
Oriented_side
side_of_plane_centered_sphereC3(const RT &ax, const RT &ay, const RT &az,
                                const RT &nx, const RT &ny, const RT &nz,
                                const RT &px, const RT &py, const RT &pz,
                                const RT &qx, const RT &qy, const RT &qz,
                                const RT &rx, const RT &ry, const RT &rz)
{
  // precondition: no two points p,q,r have the same projection
  //   <=> (p-q),(p-r), (q-r) may not be collinear to n
  //
  // resolution of the system (where we note c the center)
  //
  // ((a-c)*n)= 0             //the center lies in the plane (a, n)
  // ((c-p)* n_orthogonal)=0  // c lies in the plane containing p and n
  // (c-p)(c-p)=(c-q)(c-q)    // p,q lie on a sphere(center c)
  //
  //  return:   sign( (c-p)(c-p) - (c-r)(c-r))
  //
  //method:
  // - tranlation of p to the origin.
  // - separate computation of det and nom of the expression

  return side_of_plane_centered_sphere_translateC3<RT>(ax-px, ay-py, az-pz,
                                                       nx, ny, nz,
                                                       qx-px, qy-py,qz-pz,
                                                       rx-px, ry-py,rz-pz);
}

} //namespace CGAL

#endif // CGAL_PREDICATES_FOR_VORONOI_INTERSECTION_CARTESIAN_2_3_H
