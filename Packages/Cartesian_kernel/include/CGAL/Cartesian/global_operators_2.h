// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_GLOBAL_OPERATORS_2_H
#define CGAL_CARTESIAN_GLOBAL_OPERATORS_2_H

CGAL_BEGIN_NAMESPACE

template < class K >
inline
typename K::Point_2
operator+(const PointC2<K> &p, const VectorC2<K> &v)
{
  //return K().construct_translated_point_2_object()(p, v);
  return PointC2<K>(p.x() + v.x(), p.y() + v.y());
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
  //return K().construct_vector_2_object()(p, q);
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
  //return K().construct_scaled_vector_2_object()(w, c);
  return VectorC2<K>(w.x() * c, w.y() * c);
}

template < class K >
inline
typename K::Vector_2
operator*(const VectorC2<K> &w, const typename K::FT &c)
{
  //  return K().construct_scaled_vector_2_object()(w, c);
  return VectorC2<K>(w.x() * c, w.y() * c);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_GLOBAL_OPERATORS_2_H
