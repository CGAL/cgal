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
// $URL$
// $Id$
// 
//
// Author        : Andreas Fabri

#ifndef CGAL_CARTESIAN_VECTOR_3_H
#define CGAL_CARTESIAN_VECTOR_3_H

#include <CGAL/Origin.h>
#include <CGAL/array.h>
#include <CGAL/constant.h>

namespace CGAL {

template < class R_ >
class VectorC3
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Ray_3                Ray_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Direction_3          Direction_3;

  typedef cpp0x::array<FT, 3>               Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:

  typedef typename Rep::const_iterator      Cartesian_const_iterator;

  typedef R_                                R;

  VectorC3() {}

  VectorC3(const Null_vector &n)
  { *this = R().construct_vector_3_object()(n); }

  VectorC3(const Point_3 &a, const Point_3 &b)
  { *this = R().construct_vector_3_object()(a, b); }

  explicit VectorC3(const Segment_3 &s)
  { *this = R().construct_vector_3_object()(s); }

  explicit VectorC3(const Ray_3 &r)
  { *this = R().construct_vector_3_object()(r); }

  explicit VectorC3(const Line_3 &l)
  { *this = R().construct_vector_3_object()(l); }

  VectorC3(const FT &x, const FT &y, const FT &z)
    : base(CGAL::make_array(x, y, z)) {}

  VectorC3(const FT &x, const FT &y, const FT &z, const FT &w)
    : base( w != FT(1) ? CGAL::make_array(x/w, y/w, z/w)
                       : CGAL::make_array(x, y, z) ) {}

  const FT & x() const
  {
      return get(base)[0];
  }
  const FT & y() const
  {
      return get(base)[1];
  }
  const FT & z() const
  {
      return get(base)[2];
  }

  const FT & hx() const
  {
      return x();
  }
  const FT & hy() const
  {
      return y();
  }
  const FT & hz() const
  {
      return z();
  }
  const FT & hw() const
  {
      return constant<FT, 1>();
  }

  Cartesian_const_iterator cartesian_begin() const
  {
    return get(base).begin();
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return get(base).end();
  }

  const FT & cartesian(int i) const;
  const FT & operator[](int i) const;
  const FT & homogeneous(int i) const;

  int dimension() const
  {
      return 3;
  }

  Vector_3 operator+(const VectorC3 &w) const;
  Vector_3 operator-(const VectorC3 &w) const;
  Vector_3 operator-() const;
  Vector_3 operator/(const FT &c) const;
  FT squared_length() const;
  Direction_3 direction() const;
};

template < class R >
inline
bool
operator==(const VectorC3<R> &v, const VectorC3<R> &w)
{
  return w.x() == v.x() && w.y() == v.y() && w.z() == v.z();
}

template < class R >
inline
bool
operator!=(const VectorC3<R> &v, const VectorC3<R> &w)
{
  return !(v == w);
}

template < class R >
inline
bool
operator==(const VectorC3<R> &v, const Null_vector &) 
{
  return CGAL_NTS is_zero(v.x()) && CGAL_NTS is_zero(v.y()) &&
         CGAL_NTS is_zero(v.z());
}

template < class R >
inline
bool
operator==(const Null_vector &n, const VectorC3<R> &v) 
{
  return v == n;
}

template < class R >
inline
bool
operator!=(const VectorC3<R> &v, const Null_vector &n)
{
  return !(v == n);
}

template < class R >
inline
bool
operator!=(const Null_vector &n, const VectorC3<R> &v)
{
  return !(v == n);
}

template < class R >
inline
const typename VectorC3<R>::FT &
VectorC3<R>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i>=0) & (i<3) );
  if (i==0) return x();
  if (i==1) return y();
  return z();
}

template < class R >
inline
const typename VectorC3<R>::FT &
VectorC3<R>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
const typename VectorC3<R>::FT &
VectorC3<R>::homogeneous(int i) const
{
  if (i==3) return hw();
  return cartesian(i);
}

template < class R >
inline
typename VectorC3<R>::Vector_3
VectorC3<R>::
operator+(const VectorC3<R> &w) const
{
  return VectorC3<R>(x() + w.x(), y() + w.y(), z() + w.z());
}

template < class R >
inline
typename VectorC3<R>::Vector_3
VectorC3<R>::operator-(const VectorC3<R> &w) const
{
  return VectorC3<R>(x() - w.x(), y() - w.y(), z() - w.z());
}

template < class R >
inline
typename VectorC3<R>::Vector_3
VectorC3<R>::operator-() const
{
  return R().construct_opposite_vector_3_object()(*this);
}

template < class R >
inline
typename VectorC3<R>::FT
VectorC3<R>::squared_length() const
{
  return CGAL_NTS square(x()) + CGAL_NTS square(y()) + CGAL_NTS square(z());
}

template < class R >
inline
typename VectorC3<R>::Vector_3
VectorC3<R>::
operator/(const typename VectorC3<R>::FT &c) const
{
  return VectorC3<R>(x()/c, y()/c, z()/c);
}

template < class R >
inline
typename VectorC3<R>::Direction_3
VectorC3<R>::direction() const
{
  return Direction_3(*this);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_VECTOR_3_H
