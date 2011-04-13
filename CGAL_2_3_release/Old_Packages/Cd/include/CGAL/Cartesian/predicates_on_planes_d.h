// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/predicates_on_planes_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_PREDICATES_ON_PLANES_D_H
#define CGAL_CARTESIAN_PREDICATES_ON_PLANES_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/predicates/kernel_ftCd.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Oriented_side
side_of_oriented_plane(const PlaneCd<R CGAL_CTAG> &h,
                       const PointCd<R CGAL_CTAG> &p)
{ 
  CGAL_kernel_precondition( h.dimension() == p.dimension() );
  return side_of_oriented_planeCd(h.begin(),h.end(),p.begin());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_PLANES_3_H
