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
// file          : include/CGAL/Cartesian/global_operators_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_GLOBAL_OPERATORS_3_H
#define CGAL_CARTESIAN_GLOBAL_OPERATORS_3_H

#include <CGAL/Cartesian/redefine_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
PointC3<R CGAL_CTAG>
operator+(const PointC3<R CGAL_CTAG> &p, const VectorC3<R CGAL_CTAG> &v)
{ // FIXME : construction
  return PointC3<R CGAL_CTAG>(p.x() + v.x(), p.y() + v.y(), p.z() + v.z());
}

template < class R >
inline
PointC3<R CGAL_CTAG>
operator-(const PointC3<R CGAL_CTAG> &p, const VectorC3<R CGAL_CTAG> &v)
{ // FIXME : construction
  return PointC3<R CGAL_CTAG>(p.x() - v.x(), p.y() - v.y(), p.z() - v.z());
}

template < class R >
inline
PointC3<R CGAL_CTAG>
operator+(const Origin &, const VectorC3<R CGAL_CTAG> &v)
{
  return PointC3<R CGAL_CTAG>(v);
}

template < class R >
inline
PointC3<R CGAL_CTAG>
operator-(const Origin &, const VectorC3<R CGAL_CTAG> &v)
{
  return PointC3<R CGAL_CTAG>(-v);
}

template < class R >
inline
VectorC3<R CGAL_CTAG>
operator-(const PointC3<R CGAL_CTAG> &p, const PointC3<R CGAL_CTAG> &q)
{ // FIXME : construction
  return VectorC3<R CGAL_CTAG>(p.x() - q.x(), p.y() - q.y(), p.z() - q.z());
}

template < class R >
inline
VectorC3<R CGAL_CTAG>
operator-(const PointC3<R CGAL_CTAG> &p, const Origin &)
{
  return VectorC3<R CGAL_CTAG>(p);
}

template < class R >
inline
VectorC3<R CGAL_CTAG>
operator-(const Origin &, const PointC3<R CGAL_CTAG> &p)
{ // FIXME : construction
  return VectorC3<R CGAL_CTAG>(-p.x(), -p.y(), -p.z());
}

template < class R >
CGAL_KERNEL_INLINE
VectorC3<R CGAL_CTAG>
operator*(const typename R::FT &c, const VectorC3<R CGAL_CTAG> &w)
{ // FIXME : construction
   return VectorC3<R CGAL_CTAG>(c * w.x(), c * w.y(), c * w.z());
}

template < class R >
CGAL_KERNEL_INLINE
VectorC3<R CGAL_CTAG>
operator*(const VectorC3<R CGAL_CTAG> &w, const typename R::FT &c)
{ // FIXME : construction
   return VectorC3<R CGAL_CTAG>(c * w.x(), c * w.y(), c * w.z());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_GLOBAL_OPERATORS_3_H
