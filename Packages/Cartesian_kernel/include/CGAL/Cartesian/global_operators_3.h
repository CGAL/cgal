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

CGAL_BEGIN_NAMESPACE

template < class R >
inline
typename R::Point_3
operator+(const PointC3<R> &p, const VectorC3<R> &v)
{
  return PointC3<R>(p.x() + v.x(), p.y() + v.y(), p.z() + v.z());
}

template < class R >
inline
typename R::Point_3
operator-(const PointC3<R> &p, const VectorC3<R> &v)
{
  return PointC3<R>(p.x() - v.x(), p.y() - v.y(), p.z() - v.z());
}

template < class R >
inline
typename R::Point_3
operator+(const Origin &, const VectorC3<R> &v)
{
  return PointC3<R>(v);
}

template < class R >
inline
typename R::Point_3
operator-(const Origin &, const VectorC3<R> &v)
{
  return PointC3<R>(-v);
}

template < class R >
inline
typename R::Vector_3
operator-(const PointC3<R> &p, const PointC3<R> &q)
{
  return VectorC3<R>(p.x() - q.x(), p.y() - q.y(), p.z() - q.z());
}

template < class R >
inline
typename R::Vector_3
operator-(const PointC3<R> &p, const Origin &)
{
  return VectorC3<R>(p);
}

template < class R >
inline
typename R::Vector_3
operator-(const Origin &, const PointC3<R> &p)
{
  return VectorC3<R>(-p.x(), -p.y(), -p.z());
}

template < class R >
CGAL_KERNEL_INLINE
typename R::Vector_3
operator*(const typename R::FT &c, const VectorC3<R> &w)
{
   return VectorC3<R>(c * w.x(), c * w.y(), c * w.z());
}

template < class R >
CGAL_KERNEL_INLINE
typename R::Vector_3
operator*(const VectorC3<R> &w, const typename R::FT &c)
{
   return VectorC3<R>(c * w.x(), c * w.y(), c * w.z());
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_GLOBAL_OPERATORS_3_H
