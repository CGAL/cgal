// ==========================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/constructions_on_lines_3.h
// source        : include/CGAL/Cartesian/constructions_on_lines_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================


#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_3_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

#ifndef CGAL_CARTESIAN_POINT_3_H
#include <CGAL/Cartesian/Point_3.h>
#endif // CGAL_CARTESIAN_POINT_3_H
#ifndef CGAL_CONSTRUCTIONS_KERNEL_FTC3_H
#include <CGAL/constructions/kernel_ftC3.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
point_on_line(const LineC3<R CGAL_CTAG>& l, int i)
{
  typename R::FT x, y, z;
  point_on_lineC3(l.point().x(),l.point().y(),l.point().z(),
      l.second_point().x(),l.second_point().y(),l.second_point().z(),
      i,x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointC3<R CGAL_CTAG>
projection_line(const PointC3<R CGAL_CTAG>& p,
                const LineC3<R CGAL_CTAG>& l)
{
  typename R::FT x,y,z;
  projection_lineC3(l.point().x(), l.point().y(), l.point().z(),
      l.second_point().x(), l.second_point().y(), l.second_point().z(),
      p.x(),p.y(),p.z(),
      x,y,z);
  return PointC3<R CGAL_CTAG>(x,y,z);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_3_H
