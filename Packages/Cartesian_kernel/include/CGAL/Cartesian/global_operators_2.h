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
// file          : include/CGAL/Cartesian/global_operators_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_GLOBAL_OPERATORS_2_H
#define CGAL_CARTESIAN_GLOBAL_OPERATORS_2_H

CGAL_BEGIN_NAMESPACE

template < class K >
inline
typename K::Point_2
operator+(const PointC2<K> &p, const VectorC2<K> &v)
{
  return K().construct_translated_point_2_object()(p, v);
}

template < class K >
inline
typename K::Point_2
operator-(const PointC2<K> &p, const VectorC2<K> &v)
{
  return PointC2<K>(p.x() - v.x(), p.y() - v.y());
}

template < class K >
inline
typename K::Point_2
operator+(const Origin &, const VectorC2<K> &v)
{
  return PointC2<K>(v);
}

template < class K >
inline
typename K::Point_2
operator-(const Origin &, const VectorC2<K> &v)
{
  return PointC2<K>(-v);
}

template < class K >
inline
typename K::Vector_2
operator-(const PointC2<K> &p, const PointC2<K> &q)
{
  return VectorC2<K>(p.x() - q.x(), p.y() - q.y());
}

template < class K >
inline
typename K::Vector_2
operator-(const PointC2<K> &p, const Origin &)
{
  return VectorC2<K>(p);
}

template < class K >
inline
typename K::Vector_2
operator-(const Origin &, const PointC2<K> &p)
{
  return VectorC2<K>(-p.x(), -p.y());
}

template < class K >
inline
typename K::Vector_2
operator*(const typename K::FT &c, const VectorC2<K> &w)
{
  return K().construct_scaled_vector_2_object()(w, c);
}

template < class K >
inline
typename K::Vector_2
operator*(const VectorC2<K> &w, const typename K::FT &c)
{
  return K().construct_scaled_vector_2_object()(w, c);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_GLOBAL_OPERATORS_2_H
