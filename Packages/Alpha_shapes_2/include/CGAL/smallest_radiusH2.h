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
// file          : include/CGAL/smallest_radiusH2.h
// package       : Alpha_shapes_2(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_SMALLEST_RADIUSH2_H
#define CGAL_SMALLEST_RADIUSH2_H

#include <CGAL/PVDH2.h>

template < class FT, class RT>
FT
inline squared_radius_smallest_circumcircle(const PointH2<FT,RT> &p,
					    const PointH2<FT,RT> &q,
					    const PointH2<FT,RT> &r) const
{
  // compute the smallest radius directly
  if (orientation(p, q, r) == COLLINEAR)
    // what do we do 
    return R_FT_return(R)(0);
  else
    
    {
      Circle_2<R> c(p, q, r);
      return R_FT_return(R)(c.squared_radius());
    }
}

#endif //CGAL_SMALLEST_RADIUSH2_H
