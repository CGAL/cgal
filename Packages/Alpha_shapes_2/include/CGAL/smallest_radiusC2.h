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
#include <CGAL/PointC2.h>
#include <CGAL/VectorC2.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class FT >
inline
FT
squared_radius_smallest_circumcircle(const PointC2<FT> &p,
				     const PointC2<FT> &q,
				     const PointC2<FT> &r)
{
  
  // the computation of the squared radius takes 17 multiplications
  // and 12 additions

  VectorC2<FT> v(p-q);
  FT numerator = v*v;
  v=q-r;
  numerator *= (v*v);
  v=r-p;
  numerator *= (v*v);
 
  FT denominator = (p.x()-q.x())*(q.y()-r.y())
    - (p.y()-q.y())*(q.x()-r.x());
  
  // assert(demominator ==  (p.x()*q.y() + q.x()*r.y() + r.x()*p.y()
  //			  -r.x()*q.y() - q.x()*p.y() - p.x()*r.y()));

  return ((denominator > FT(0) ? 
	   numerator /(FT(4) * denominator * denominator) : FT(0)));
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_SMALLEST_RADIUSC2_H
