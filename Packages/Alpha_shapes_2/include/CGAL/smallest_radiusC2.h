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
// file          : include/CGAL/smallest_radiusC2.h
// package       : Alpha_shapes_2(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
//
// ======================================================================

#ifndef CGAL_SMALLEST_RADIUSC2_H
#define CGAL_SMALLEST_RADIUSC2_H

#include <CGAL/cartesian_classes.h>
//#include <CGAL/Cartesian/Point_2.h>
//#include <CGAL/Cartesian/Vector_2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class R >
inline
typename R::FT
squared_radius_smallest_circumcircle(const PointC2<R> &p,
				     const PointC2<R> &q,
				     const PointC2<R> &r)
{
  
  // the computation of the squared radius takes 17 multiplications
  // and 12 additions
  
  typedef typename R::FT FTT;

  VectorC2<R> v(p-q);
  FTT numerator = v*v;
  v=q-r;
  numerator *= (v*v);
  v=r-p;
  numerator *= (v*v);
 
  FTT denominator = (p.x()-q.x())*(q.y()-r.y())
    - (p.y()-q.y())*(q.x()-r.x());
  
  // assert(demominator ==  (p.x()*q.y() + q.x()*r.y() + r.x()*p.y()
  //			  -r.x()*q.y() - q.x()*p.y() - p.x()*r.y()));

  return ((denominator > FTT(0) ? 
	   numerator /(FTT(4) * denominator * denominator) : FTT(0)));
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_SMALLEST_RADIUSC2_H
