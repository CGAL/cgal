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
// file          : include/CGAL/Cartesian/constructions_on_lines_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_D_H
#define CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/Cartesian/Point_d.h>
#include <CGAL/Cartesian/Line_d.h>
#include <CGAL/constructions/kernel_ftCd.h>

CGAL_BEGIN_NAMESPACE

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointCd<R CGAL_CTAG>
point_on_line(typename R::FT ratio, const LineCd<R CGAL_CTAG>& l)
{
  PointCd<R CGAL_CTAG> p(l.dimension());
  const PointCd<R CGAL_CTAG> lp(l.point());
  const DirectionCd<R CGAL_CTAG> ld(l.direction());
  point_on_lineCd(lp.begin(), lp.end(), ld.begin(), ratio,
		  p.begin(), R());
  return p;
}

template <class R>
CGAL_KERNEL_LARGE_INLINE
PointCd<R CGAL_CTAG>
projection_line(const PointCd<R CGAL_CTAG>& p,
                const LineCd<R CGAL_CTAG>& l)
{
  PointCd<R CGAL_CTAG> q(p.dimension());
  const PointCd<R CGAL_CTAG> lp(l.point());
  const DirectionCd<R CGAL_CTAG> ld(l.direction());
  projection_lineCd(p.begin(),p.end(),
		    lp.begin(), lp.end(), ld.begin(),ld.end(),
                    q.begin(), R());
  return q;
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_D_H
