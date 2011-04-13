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
// file          : include/CGAL/Cartesian/predicates_on_planes_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_PREDICATES_ON_PLANES_3_H
#define CGAL_CARTESIAN_PREDICATES_ON_PLANES_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/predicates/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Oriented_side
side_of_oriented_plane(const PlaneC3<R CGAL_CTAG> &h,
                       const PointC3<R CGAL_CTAG> &p)
{ 
  return side_of_oriented_planeC3(h.a(), h.b(), h.c(), h.d(),
	                          p.x(), p.y(), p.z());
}

template < class R >
inline
bool
equal_plane(const PlaneC3<R CGAL_CTAG> &h, const PlaneC3<R CGAL_CTAG> &p)
{ 
  return equal_planeC3(h.a(), h.b(), h.c(), h.d(),
	               p.a(), p.b(), p.c(), p.d());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_PLANES_3_H
