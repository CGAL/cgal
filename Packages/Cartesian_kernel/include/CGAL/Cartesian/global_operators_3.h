// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri, Herve Bronnimann

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
  return PointC3<K>(v.x(), v.y(), v.z());
}

template < class K >
inline
typename K::Point_3
operator-(const Origin &, const VectorC3<K> &v)
{
  return PointC3<K>(-v.x(), -v.y(), -v.z());
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
  return VectorC3<K>(p.x(), p.y(), p.z());
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
