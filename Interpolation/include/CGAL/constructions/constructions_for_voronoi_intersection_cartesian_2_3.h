// Copyright (c) 2003 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Julia Floetotto

#ifndef CGAL_CONSTRUCTIONS_FOR_VORONOI_INTERSECTION_CARTESIAN_2_3_H
#define CGAL_CONSTRUCTIONS_FOR_VORONOI_INTERSECTION_CARTESIAN_2_3_H

#include <CGAL/license/Interpolation.h>


namespace CGAL {

template < class RT>
void
plane_centered_circumcenter_translateC3(const RT &ax, const RT &ay,
                                        const RT &az,
                                        const RT &nx, const RT &ny,
                                        const RT &nz,
                                        const RT &qx, const RT &qy,
                                        const RT &qz,
                                        const RT &rx, const RT &ry,
                                        const RT &rz,
                                        RT &x, RT &y, RT &z)
{
  RT den = RT(2) * determinant(nx,qx,rx,
                               ny,qy,ry,
                               nz,qz,rz);
  // The 3 points aren't collinear.
  // Hopefully, this is already checked at the upper level.
  CGAL_assertion ( den != RT(0) );

  RT q2 = CGAL_NTS square(qx) + CGAL_NTS square(qy) +
          CGAL_NTS square(qz);
  RT r2 = CGAL_NTS square(rx) + CGAL_NTS square(ry) +
          CGAL_NTS square(rz);
  RT na = nx*ax + ny*ay + nz*az;
  na *= RT(2.0);

  x =   determinant(ny,nz,na,
                    qy,qz,q2,
                    ry,rz,r2)/ den ;

  y = - determinant(nx,nz,na,
                    qx,qz,q2,
                    rx,rz,r2)/ den ;

  z =   determinant(nx,ny,na,
                    qx,qy,q2,
                    rx,ry,r2)/ den ;
}

template < class RT>
void
plane_centered_circumcenterC3(const RT &ax, const RT &ay, const RT &az,
                              const RT &nx, const RT &ny, const RT &nz,
                              const RT &px, const RT &py, const RT &pz,
                              const RT &qx, const RT &qy, const RT &qz,
                              const RT &rx, const RT &ry, const RT &rz,
                              RT &x, RT &y, RT &z)
{
  // resolution of the system (where we note c the center)
  //
  // ((a-c)*n)= 0                //the center lies in the plane (a, n)
  // (c-p)(c-p)=(c-q)(c-q)=(c-r)(c-r)   // p,q,r on a sphere(center c)
  //
  //precondition: p,q,r aren't collinear.
  //method:
  // - tranlation of p to the origin.
  plane_centered_circumcenter_translateC3(ax-px, ay-py, az-pz,
                                          nx, ny, nz,
                                          qx-px, qy-py,qz-pz,
                                          rx-px, ry-py,rz-pz,
                                          x, y, z);
  x+=px;
  y+=py;
  z+=pz;
}

template < class RT>
void
bisector_plane_intersection_translateC3(const RT &ax, const RT &ay,
                                        const RT &az,
                                        const RT &nx, const RT &ny,
                                        const RT &nz,
                                        const RT &qx, const RT &qy,
                                        const RT &qz, const RT& den,
                                        RT &x1, RT &y1, RT &x2, RT
                                        &y2, bool& swapped)
{
  // c: a point on l must be the center of a sphere passing
  // through p and q, c lies in h. 2 equations:
  //   c^2 = (q-c)^2   //p and q on the sphere's boundary
  //  (c-a)n = 0       //c in the plane h
  // the line is defined by p1 and p2 with
  //=> p1: z1 =0, p2: z2=1
  // precondition: (nx!=0 || ny!=0) && (qx!=0 && qy!=0) && den!=0
  // where RT den = RT(2.0) * determinant(qx,qy,nx, ny);

  RT q2 = CGAL_NTS square(qx) + CGAL_NTS square(qy)
          + CGAL_NTS square(qz);
  RT na = nx*ax + ny*ay + nz*az;
  na *= RT(2.0);

  x1 = determinant(ny, na, qy, q2);
  y1 = - determinant(nx, na, qx, q2);

  x2 = x1 +  RT(2.0) * determinant(qy,qz,ny, nz);
  y2 = y1 -  RT(2.0) * determinant(qx,qz,nx, nz);

  x1 /= den;
  x2 /= den;
  y1 /= den;
  y2 /= den;

  //we need to orient the line such that
  // if p is on the positive side of the plane
  //(<=> (p-a)*n >0 <=> na < 0) then orientation (pq p1 p2) is ccw
  // if not: permutation of p1 and p2
  if((sign_of_determinant(qx,qy,qz, x1,y1,RT(0),x2 ,y2,RT(1))
      * CGAL_NTS sign (-na)) > 0 )
  {
    RT x3(x1),y3(y1);
    x1 =x2;
    y1 =y2;
    x2 = x3;
    y2 = y3;
    swapped =true;
  }
}

template < class RT>
void
bisector_plane_intersection_permuteC3(const RT &ax, const RT &ay,
                                      const RT &az,
                                      const RT &nx, const RT &ny,
                                      const RT &nz,
                                      const RT &px, const RT &py,
                                      const RT &pz,
                                      const RT &qx, const RT &qy,
                                      const RT &qz,
                                      const RT &den,
                                      RT &x1, RT &y1, RT& z1,
                                      RT &x2, RT &y2, RT& z2)
{
  //translation of p to the origin
  bool swapped =false;
  CGAL_precondition((nx!=RT(0) ||  ny!=RT(0)) && (qx!=px || qy!=py)
                    &&den!=RT(0));

  bisector_plane_intersection_translateC3(ax-px, ay-py, az-pz,
                                          nx, ny, nz,
                                          qx-px, qy-py,qz-pz,den,
                                          x1, y1,x2,y2,swapped);
  // re-translation of the origin to p:
  x1+=px;
  y1+=py;
  z1 = pz;
  x2+=px;
  y2+=py;
  z2 = pz+RT(1);
  if (swapped){
    z1 += RT(1);
    z2 =pz;
  }
}

template < class RT>
void
bisector_plane_intersectionC3(const RT &ax, const RT &ay, const RT &az,
                              const RT &nx, const RT &ny, const RT &nz,
                              const RT &px, const RT &py, const RT &pz,
                              const RT &qx, const RT &qy, const RT &qz,
                              RT &x1, RT &y1, RT& z1,
                              RT &x2, RT &y2, RT& z2)
{
  // constructs the line l = (p1,p2)= ((x1,y1,z1),(x2,y2,z2))
  // the intersection line between the bisector of (p,q) and
  //  h: plane containing a orthogonal to n.
  // l is oriented such that if p is on the positive side of h,
  // then (pq p1 p2) is ccw, otherwise (pq p1 p2) is cw
  //
  // precondition: (pq) is not parallel to n
  //
  // method: computing the intersection points of l
  //    with the planes (z=0/z=1) or (x=0/x=1) or (y=0/y=1) depending if
  //    1) the line is not parallel to the plane
  //    2) the projection of n and (p-q) onto the plane is not
  //    identical
  //    computation for (z=0) with adequate permutations
  RT den = RT(2.0) * determinant(qx-px,qy-py,nx, ny);
  if ((nx!=0 ||  ny!=0) && (qx!=px || qy!=py) && den!=RT(0))
    //den==0 <=> projections of (qx,qy) and (nx,ny) are identical
    //intersection with z=0/z=1
    bisector_plane_intersection_permuteC3(ax,ay,az,nx,ny,nz,px,py,pz,
                                          qx,qy,qz,den,
                                          x1,y1,z1,x2,y2,z2);
  else{
    den = RT(2.0) * determinant(qy-py,qz-pz,ny,nz);
    if ((ny!=0 ||  nz!=0) && (qy!=py || qz!=pz) && den!=RT(0))
      //intersection with x=0/x=1 => permutations
      bisector_plane_intersection_permuteC3(ay,az,ax,ny,nz,nx,py,pz,px,
                                            qy,qz,qx,den,
                                            y1,z1,x1,y2,z2,x2);
    else{
      den = RT(2.0) * determinant(qz-pz,qx-px,nz,nx);
      CGAL_assertion((nx!=0 ||  nz!=0) && (qx!=px || qz!=pz) && den!=RT(0));
      //intersection with y=0/y=1 => permutations
      bisector_plane_intersection_permuteC3(az,ax,ay,nz,nx,ny,pz,px,py,
                                            qz,qx,qy,den,
                                            z1,x1,y1,z2,x2,y2);
    }
  }
}

} //namespace CGAL

#endif // CGAL_CONSTRUCTIONS_FOR_VORONOI_INTERSECTION_CARTESIAN_2_3_H
