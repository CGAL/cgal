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
// file          : include/CGAL/smallest_radius_ftC3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef SMALLEST_RADIUSFT3_H
#define SMALLEST_RADIUSFT3_H

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class FT >
FT
squared_radius_smallest_circumsphere(const FT &px, const FT &py, const FT &pz,
				     const FT &qx, const FT &qy, const FT &qz)
{
  FT vx = px-qx;
  FT vy = py-qy;
  FT vz = pz-qz;
  return (vx*vx + vy*vy + vz*vz)/FT(4);
}

//---------------------------------------------------------------------

template < class FT >
FT
squared_radius_smallest_circumsphere(const FT &px, const FT &py, const FT &pz,
				     const FT &qx, const FT &qy, const FT &qz,
				     const FT &rx, const FT &ry, const FT &rz)
{
  // the computation of the squared radius takes 22 multiplications
  // and x additions

  FT vx = px-qx;
  FT vy = py-qy;
  FT vz = pz-qz;
  FT numerator = vx*vx + vy*vy + vz*vz;
  vx=qx-rx;
  vy=qy-ry;
  vz=qz-rz;
  numerator *= (vx*vx + vy*vy + vz*vz);
  vx=rx-px;
  vy=ry-py;
  vz=rz-pz;
  numerator *= (vx*vx + vy*vy + vz*vz);

  FT det = (px-qx)*(qy-ry)
         - (py-qy)*(qx-rx);


  FT denominator = det*det;
  
  det = (px-qx)*(qz-rz)
    - (pz-qz)*(qx-rx);

  denominator += (det*det);

  det = (py-qy)*(qz-rz)
    - (pz-qz)*(qy-ry);

  denominator += (det*det);
 
  return ((denominator > FT0 ? 
		       numerator /(FT(4) * denominator) : FT0));
}

//---------------------------------------------------------------------

template < class FT >
FT
squared_radius_circumsphere(const FT &px, const FT &py, const FT &pz, 
			    const FT &qx, const FT &qy, const FT &qz,
			    const FT &rx, const FT &ry, const FT &rz,
			    const FT &sx, const FT &sy, const FT &sz)
{
  FT FT0(0);
  FT FT4(4);

  FT p3 = px*px + py*py + pz*pz; 
  FT q3 = qx*qx + qy*qy + qz*qz; 
  FT r3 = rx*rx + ry*ry + rz*rz;
  FT s3 = sx*sx + sy*sy + sz*sz;

  // Maple's result takes only 84 multiplications

  FT t1,t10,t11,t12,t13,t15,t16,t17,t19,t2,t22;
  FT t28,t29,t3,t30,t31,t32,t33,t34,t35,t37,t38;
  FT t39,t4,t40,t41,t43,t44,t45,t47,t5,t50;
  FT t6,t68,t7,t81,t9;

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
  FT det_sphere = pz*t35 - qz*t41 + rz*t45 - sz*t47;
  if (det_sphere == FT0)
    return FT0;
  FT det = p3*t35 - q3*t41 + r3*t45 - s3*t47;
  t50 = det*det;

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

  det = py*t7 - qy*t13 + ry*t17 - sy*t19;
  t22 = det*det;
  det = px*t7 - qx*t13 + rx*t17 - sx*t19;
  t28 = det*det;
  
  t68 = t43*t4-t43*t6-t37*t2+t37*t5+t38*t1-t38*t3-t44*t4+t44*t6+t29*t10-t29
    *t12-t30*t9+t30*t11;
  t81 = t39*t2-t39*t5-t31*t10+t31*t12+t32*t15-t32*t16-t40*t1+t40*t3+t33*t9-
    t33*t11-t34*t15+t34*t16;

  
  FT result = (t22+t28+t50+FT4*det_sphere*(t68+t81)) / 
              (FT4 * det_sphere * det_sphere);

  return result;
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //SMALLEST_RADIUSFT3_H
