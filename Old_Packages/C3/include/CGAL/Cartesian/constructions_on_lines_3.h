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
// file          : include/CGAL/Cartesian/constructions_on_lines_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_3_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/constructions/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
point_on_line(int i, const LineC3<R CGAL_CTAG>& l)
{
  typename R::FT x, y, z;
  point_on_lineC3(l.point().x(),l.point().y(),l.point().z(),
                  l.direction().dx(),l.direction().dy(),l.direction().dz(),
                  i,x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
projection_line(const PointC3<R CGAL_CTAG>& p, const LineC3<R CGAL_CTAG>& l)
{
  typename R::FT x,y,z;
  projection_lineC3(p.x(),p.y(),p.z(),
		    l.point().x(), l.point().y(), l.point().z(),
                    l.direction().dx(),l.direction().dy(),l.direction().dz(),
                    x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_3_H
