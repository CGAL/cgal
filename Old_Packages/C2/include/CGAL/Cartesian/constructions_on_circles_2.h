// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/constructions_on_circles_2.h
// source        : include/CGAL/Cartesian/constructions_on_circles_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann (hbronni@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_CIRCLES2_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_CIRCLES2_H

#include <utility>

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#include <CGAL/Cartesian/redefine_names_2.h>
#endif

#ifndef CGAL_CARTESIAN_POINT_2_H
#include <CGAL/Cartesian/Point_2.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
inline
typename R::FT
squared_circumradius( PointC2<R CGAL_CTAG> const& p,
                      PointC2<R CGAL_CTAG> const& q,
                      PointC2<R CGAL_CTAG> const& r)
{
  return squared_circumradiusC2(p.x(),p.y(),q.x(),q.y(),r.x(),r.y());
}

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

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_CIRCLES_2_H
