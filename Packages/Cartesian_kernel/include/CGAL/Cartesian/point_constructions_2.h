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
// file          : include/CGAL/Cartesian/point_constructions_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_POINT_CONSTRUCTIONS_2_H
#define CGAL_CARTESIAN_POINT_CONSTRUCTIONS_2_H

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

template < class R >
inline
PointC2<R>
line_get_point( LineC2<R CGAL_CTAG> const& l, int i)
{
  typename R::FT x, y;
  line_get_pointC2(l.a(), l.b(), l.c(), i, x, y);
  return PointC2<R>(x,y);
}

template < class R >
inline
PointC2<R>
line_project_point( LineC2<R CGAL_CTAG> const& l,
                    PointC2<R CGAL_CTAG> const& p)
{
  typename R::FT x, y;
  line_project_pointC2(l.a(), l.b(), l.c(), p.x(), p.y(), x, y);
  return PointC2<R>(x,y);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_CONSTRUCTIONS_2_H
