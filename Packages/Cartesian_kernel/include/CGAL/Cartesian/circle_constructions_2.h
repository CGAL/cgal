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
// file          : include/CGAL/Cartesian/circle_constructions_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_CIRCLE_CONSTRUCTIONS_2_H
#define CGAL_CARTESIAN_CIRCLE_CONSTRUCTIONS_2_H

#include <utility>
#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/Cartesian/Point_2.h>

CGAL_BEGIN_NAMESPACE

// Ok this function should not be here.
template < class R >
inline
typename R::FT
squared_circumradius( PointC2<R CGAL_CTAG> const& p,
                      PointC2<R CGAL_CTAG> const& q,
                      PointC2<R CGAL_CTAG> const& r)
{
  return squared_circumradiusC2(p.x(),p.y(),q.x(),q.y(),r.x(),r.y());
}

// This one's used nowhere.  Probably should return a Circle ?
template < class R >
inline
PointC2<R CGAL_CTAG>
squared_circumcircle( PointC2<R CGAL_CTAG> const& p,
                      PointC2<R CGAL_CTAG> const& q,
                      PointC2<R CGAL_CTAG> const& r,
                      typename R::FT &radius)
{
  typename R::FT x,y;
  radius = squared_circumradiusC2(p.x(),p.y(),q.x(),q.y(),r.x(),r.y(),x,y);
  return PointC2<R CGAL_CTAG>(x,y);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CIRCLE_CONSTRUCTIONS_2_H
