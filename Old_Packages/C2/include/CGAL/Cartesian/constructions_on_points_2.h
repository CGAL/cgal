// ============================================================================
//
// Copyright (c) 1998, 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/constructions_on_points_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================

#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_POINTS_2_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_POINTS_2_H

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/Cartesian/Point_2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
PointC2<R CGAL_CTAG>
midpoint( PointC2<R CGAL_CTAG> const& p,
          PointC2<R CGAL_CTAG> const& q )
{
  typename R::FT x,y;
  midpointC2(p.x(),p.y(),q.x(),q.y(),x,y);
  return PointC2<R CGAL_CTAG>(x,y);
}

template < class R >
inline
PointC2<R CGAL_CTAG>
circumcenter( PointC2<R CGAL_CTAG> const& p,
              PointC2<R CGAL_CTAG> const& q,
              PointC2<R CGAL_CTAG> const& r)
{
  typename R::FT x,y;
  circumcenterC2(p.x(),p.y(),q.x(),q.y(),r.x(),r.y(),x,y);
  return PointC2<R CGAL_CTAG>(x,y);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_POINTS_2_H
