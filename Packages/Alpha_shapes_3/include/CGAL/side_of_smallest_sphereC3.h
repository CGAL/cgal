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
// file          : include/CGAL/side_of_smallest_sphereC3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef SIDE_OF_SMALLEST_SPHEREC3_H
#define SIDE_OF_SMALLEST_SPHEREC3_H

#include <CGAL/cartesian_classes.h>
//#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/predicates/side_of_smallest_sphere_ftC3.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class R >
inline
Bounded_side
side_of_bounded_sphere(const PointC3<R> &p,
		       const PointC3<R> &q,
		       const PointC3<R> &r,
		       const PointC3<R> &test) 
{

  return side_of_bounded_sphere(p.x(), p.y(), p.z(),
				q.x(), q.y(), q.z(),
				r.x(), r.y(), r.z(),
				test.x(), test.y(), test.z());

  // the direct computation of the squared radius takes 213 multiplications
  // and x additions
  // developping 4x4 determinant already takes 40 multiplications
  // whereas /Maple's result takes only 64 multiplications

  /*
  FT det4 = det4x4_by_formula(p[1], p[2], p3, FT1,
				   q[1], q[2], q3, FT1,
				   r[1], r[2], r3, FT1,
				   s[1], s[2], s3, FT1);

  FT det3 = det3x3_by_formula(p[1], p[2], FT1,
				   q[1], q[2], FT1,
				   r[1], r[2], FT1);
	
  FT result1 = det4 * det3;

  det4 = det4x4_by_formula(p[0], p[2], p3, FT1,
				q[0], q[2], q3, FT1,
				r[0], r[2], r3, FT1,
				s[0], s[2], s3, FT1);

  det3 = det3x3_by_formula(p[0], p[2], FT1,
				q[0], q[2], FT1,
				r[0], r[2], FT1);

  result1 += det4 * det3;

  det4 = det4x4_by_formula(p[0], p[1], p3, FT1,
				q[0], q[1], q3, FT1,
				r[0], r[1], r3, FT1,
				s[0], s[1], s3, FT1);
  
  det3 = det3x3_by_formula(p[0], p[1], FT1,
				q[0], q[1], FT1,
				r[0], r[1], FT1);

  result1 += det4 * det3;

  det4 = det4x4_by_formula(p[0], p[1], p[2], FT1,
				q[0], q[1], q[2], FT1,
				r[0], r[1], r[2], FT1,
				s[0], s[1], s[2], FT1);

  det3 = det3x3_by_formula(p[0], p[1], p[2],
				q[0], q[1], q[2],
				r[0], r[1], r[2]);
  
  result1 -= FT(2) * det4 * det3;			  
  */
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //SIDE_OF_SMALLEST_SPHEREC3_H
