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

template < class K >
inline
typename K::Point_3
operator+(const PointC3<K> &p, const VectorC3<K> &v)
{
  return K().construct_translated_point_3_object()(p, v);
}

template < class K >
inline
typename K::Point_3
operator-(const PointC3<K> &p, const VectorC3<K> &v)
{
  return PointC3<K>(p.x() - v.x(), p.y() - v.y(), p.z() - v.z());
}

template < class K >
inline
typename K::Point_3
operator+(const Origin &, const VectorC3<K> &v)
{
  return PointC3<K>(v);
}

template < class K >
inline
typename K::Point_3
operator-(const Origin &, const VectorC3<K> &v)
{
  return PointC3<K>(-v);
}

template < class K >
inline
typename K::Vector_3
operator-(const PointC3<K> &p, const PointC3<K> &q)
{
  return VectorC3<K>(p.x() - q.x(), p.y() - q.y(), p.z() - q.z());
  //return K().construct_vector_3_object()(q, p);
}

template < class K >
inline
typename K::Vector_3
operator-(const PointC3<K> &p, const Origin &)
{
  return VectorC3<K>(p);
}

template < class K >
inline
typename K::Vector_3
operator-(const Origin &, const PointC3<K> &p)
{
  return VectorC3<K>(-p.x(), -p.y(), -p.z());
}

template < class K >
inline
typename K::Vector_3
operator*(const typename K::FT &c, const VectorC3<K> &w)
{
  return K().construct_scaled_vector_3_object()(w, c);
}

template < class K >
inline
typename K::Vector_3
operator*(const VectorC3<K> &w, const typename K::FT &c)
{
  return K().construct_scaled_vector_3_object()(w, c);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_GLOBAL_OPERATORS_3_H
