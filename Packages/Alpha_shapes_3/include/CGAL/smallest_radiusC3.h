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
// file          : include/CGAL/smallest_radiusC3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef SMALLEST_RADIUSC3_H
#define SMALLEST_RADIUSC3_H

#include <CGAL/cartesian_classes.h>
//#include <CGAL/Cartesian/Point_3.h>
//#include <CGAL/Cartesian/VectorC3.h>
#include <CGAL/constructions/smallest_radius_ftC3.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class R >
inline
typename R::FT
squared_radius_smallest_circumsphere(const PointC3<R> &p,
				     const PointC3<R> &q,
				     const PointC3<R> &r) 
{
  // the computation of the squared radius takes 22 multiplications
  // and x additions
  typedef typename R::FT FT;

  FT FT0(0);
  VectorC3<R> v(p-q);
  FT numerator = v*v;
  v=q-r;
  numerator *= (v*v);
  v=r-p;
  numerator *= (v*v);

  FT det = (p[0]-q[0])*(q[1]-r[1])
         - (p[1]-q[1])*(q[0]-r[0]);

  // cerr.precision(12);
  // cerr << det << " " << (p[0]*q[1] + q[0]*r[1] + r[0]*p[1]
  // 		        -r[0]*q[1] - q[0]*p[1] - p[0]*r[1]) << endl;

  FT denominator = det*det;
  
  det = (p[0]-q[0])*(q[2]-r[2])
    - (p[2]-q[2])*(q[0]-r[0]);

  denominator += (det*det);

  det = (p[1]-q[1])*(q[2]-r[2])
    - (p[2]-q[2])*(q[1]-r[1]);

  denominator += (det*det);
 
  return ((denominator > FT0 ? 
		       numerator /(FT(4) * denominator) : FT0));
}

template < class R >
inline
typename R::FT
squared_radius_circumsphere(const PointC3<R> &p,
			    const PointC3<R> &q,
			    const PointC3<R> &r,
			    const PointC3<R> &s) 
{
  return squared_radius_circumsphere(p.x(), p.y(), p.z(), 
				     q.x(), q.y(), q.z(),
				     r.x(), r.y(), r.z(),
				     s.x(), s.y(), s.z());

  // the direct computation of the squared radius takes 220 multiplications
  // and x additions
  // developping 4x4 determinant takes 40 multiplications
  // Maple's result takes only 84 multiplications
  /*
  FT det_sphere = det4x4_by_formula(p[0], p[1], p[2], FT1,
					 q[0], q[1], q[2], FT1,
					 r[0], r[1], r[2], FT1,
					 s[0], s[1], s[2], FT1);

  if (det_sphere == FT0)
    return FT0;

   FT det = det4x4_by_formula(p[1], p[2], p3, FT1,
				  q[1], q[2], q3, FT1,
				  r[1], r[2], r3, FT1,
				  s[1], s[2], s3, FT1);
  FT numerator = (det*det);

  det = det4x4_by_formula(p[0], p[2], p3, FT1,
			       q[0], q[2], q3, FT1,
			       r[0], r[2], r3, FT1,
			       s[0], s[2], s3, FT1);
  numerator += (det*det);

  det = det4x4_by_formula(p[0], p[1], p3, FT1,
			       q[0], q[1], q3, FT1,
			       r[0], r[1], r3, FT1,
			       s[0], s[1], s3, FT1);
  numerator += (det*det);
  
  det = det4x4_by_formula(p[0], p[1], p[2], p3,
			       q[0], q[1], q[2], q3,
			       r[0], r[1], r[2], r3,
			       s[0], s[1], s[2], s3);
  

  numerator += (FT4*det_sphere*det);

  FT result1 = numerator / (FT4 * det_sphere * det_sphere);
  */
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //SMALLEST_RADIUSC3_H
