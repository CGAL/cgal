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
// file          : include/CGAL/Cartesian/predicates_on_directions_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_PREDICATES_ON_DIRECTIONS_D_H
#define CGAL_PREDICATES_ON_DIRECTIONS_D_H

#include <CGAL/predicates/kernel_ftCd.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
bool
equal_direction(const DirectionCd<R CGAL_CTAG>& d1,
                const DirectionCd<R CGAL_CTAG>& d2)
{
  if (d1.dimension() != d2.dimension()) return false;
  return equal_directionCd(d1.begin(),d1.end(),d2.begin());
}

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_ON_DIRECTIONS_D_H
