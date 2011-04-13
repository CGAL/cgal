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
// file          : include/CGAL/predicates/in_smallest_orthogonalcircle_ftC2.h
// package       : Alpha_shapes_2 (1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_IN_SMALLEST_ORTHOGONALCIRCLE_FTC2_H 
#define CGAL_IN_SMALLEST_ORTHOGONALCIRCLE_FTC2_H

#include <CGAL/determinant.h>
#include <CGAL/enum.h>

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
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
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif //CGAL_IN_SMALLEST_ORTHOGONALCIRCLEC2_H
