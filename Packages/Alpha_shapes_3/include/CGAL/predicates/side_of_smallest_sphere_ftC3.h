// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-20 $
// release_date  : $CGAL_Date: 1999/06/02 $
//
// file          : include/CGAL/predicates/side_of_smallest_sphere_ftC3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef SIDE_OF_SMALLEST_SPHERE_FTC3_H
#define SIDE_OF_SMALLEST_SPHERE_FTC3_H

#include <CGAL/constructions/smallest_radius_ftC3.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template <class FT>
inline 
Bounded_side 
side_of_bounded_sphere(const FT &px, const FT &py, const FT &pz,
		       const FT &qx, const FT &qy, const FT &qz,
		       const FT &testx, const FT &testy, const FT &testz)
{
  FT vx = testx - (px + qx)/FT(2);
  FT vy = testy - (py + qy)/FT(2);
  FT vz = testz - (pz + qz)/FT(2);
  FT squared_distance = vx*vx + vy*vy + vz*vz;

  FT squared_radius = squared_radius_smallest_circumsphere(px, py, pz,
							   qx, qy, qz);
  
  return (squared_radius > squared_distance) ?
                      ON_BOUNDED_SIDE :
                      ((squared_radius < squared_distance) ?
                                            ON_UNBOUNDED_SIDE :
                                            ON_BOUNDARY);
}

//------------------------------------------------------------------------

template < class FT >
inline
Bounded_side
side_of_bounded_sphere(const FT &px, const FT &py, const FT &pz, 
		       const FT &qx, const FT &qy, const FT &qz,
		       const FT &rx, const FT &ry, const FT &rz,
		       const FT &sx, const FT &sy, const FT &sz)
{
  FT FT0(0);

  FT p3 = px*px + py*py + pz*pz; 
  FT q3 = qx*qx + qy*qy + qz*qz; 
  FT r3 = rx*rx + ry*ry + rz*rz;
  FT s3 = sx*sx + sy*sy + sz*sz;

   // Maple's result takes only 64 multiplications

  FT t1,t10,t11,t12,t13,t15,t16,t17,t19,t2;
  FT t29,t3,t30,t31,t32,t33,t34,t35,t37,t38;
  FT t39,t4,t40,t41,t43,t44,t45,t47,t5, t6, t7, t9;

  t1 = qz*r3;

  t2 = qz*s3;
  t3 = rz*q3;
  t4 = rz*s3;
  t5 = sz*q3;
  t6 = sz*r3;
  t7 = t1-t2-t3+t4+t5-t6;
  t9 = pz*r3;
  t10 = pz*s3;
  t11 = rz*p3;
  t12 = sz*p3;
  t13 = t9-t10-t11+t4+t12-t6;
  t15 = pz*q3;
  t16 = qz*p3;
  t17 = t15-t10-t16+t2+t12-t5;
  t19 = t15-t9-t16+t1+t11-t3;

  // developping the first column

  FT det4 = py*t7 - qy*t13 + ry*t17 - sy*t19;
  FT det3 = (py-qy)*(qz-rz) - (pz-qz)*(qy-ry);
  FT result = det4 *det3;

  det4 = px*t7 - qx*t13 + rx*t17 - sx*t19;
  det3 = (px-qx)*(qz-rz) - (pz-qz)*(qx-rx);
  result += det4 *det3;

  t29 = qx*ry;
  t30 = qx*sy;
  t31 = rx*qy;
  t32 = rx*sy;
  t33 = sx*qy;
  t34 = sx*ry;
  t35 = t29-t30-t31+t32+t33-t34;
  t37 = px*ry;
  t38 = px*sy;
  t39 = rx*py;
  t40 = sx*py;
  t41 = t37-t38-t39+t32+t40-t34;
  t43 = px*qy;
  t44 = qx*py;
  t45 = t43-t38-t44+t30+t40-t33;
  t47 = t43-t37-t44+t29+t39-t31;

  // developping the third column
  det4 = p3*t35 - q3*t41 + r3*t45 - s3*t47;
  det3 = (px-qx)*(qy-ry) - (py-qy)*(qx-rx);
  result += det4 * det3;

  det4 = pz*t35 - qz*t41 + rz*t45 - sz*t47;
  det3 = pz*(t29-t31) - qz*(t37-t39) + rz*(t43-t44);
  result -= FT(2) * det4 * det3;

  return   (result > FT0 ? ON_BOUNDED_SIDE :
	              ((result < FT0) ?
		       ON_UNBOUNDED_SIDE :
		       ON_BOUNDARY));
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif  // SIDE_OF_SMALLEST_SPHERE_FTC3_H
