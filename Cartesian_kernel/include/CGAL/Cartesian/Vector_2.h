// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_VECTOR_2_H
#define CGAL_CARTESIAN_VECTOR_2_H

#include <CGAL/Origin.h>
#include <CGAL/array.h>
#include <CGAL/constant.h>
#include <CGAL/Handle_for.h>

namespace CGAL {

template < class R_ >
class VectorC2
{
  typedef VectorC2<R_>                      Self;
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Direction_2          Direction_2;

  typedef std::array<FT, 2>               Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:

  typedef typename Rep::const_iterator      Cartesian_const_iterator;

  typedef R_                                R;

  VectorC2() {}

  VectorC2(const FT &x, const FT &y)
    : base(CGAL::make_array(x, y)) {}

  VectorC2(const FT &hx, const FT &hy, const FT &hw)
    : base( hw != FT(1) ? CGAL::make_array<FT>(hx/hw, hy/hw)
                        : CGAL::make_array(hx, hy) ) {}

  friend void swap(Self& a, Self& b)
#ifdef __cpp_lib_is_swappable
    noexcept(std::is_nothrow_swappable_v<Base>)
#endif
  {
    using std::swap;
    swap(a.base, b.base);
  }

  const FT & x() const
  {
      return CGAL::get_pointee_or_identity(base)[0];
  }

  const FT & y() const
  {
      return CGAL::get_pointee_or_identity(base)[1];
  }

  const FT & hx() const
  {
      return x();
  }

  const FT & hy() const
  {
      return y();
  }

  const FT& hw() const
  {
    return constant<FT, 1>();
  }

  Cartesian_const_iterator cartesian_begin() const
  {
    return CGAL::get_pointee_or_identity(base).begin();
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return CGAL::get_pointee_or_identity(base).end();
  }

};

template < class R >
CGAL_KERNEL_INLINE
bool
operator==(const VectorC2<R> &v, const VectorC2<R> &w)
{
  return w.x() == v.x() && w.y() == v.y();
}

template < class R >
inline
bool
operator!=(const VectorC2<R> &v, const VectorC2<R> &w)
{
  return !(v == w);
}

template < class R >
inline
bool
operator==(const VectorC2<R> &v, const Null_vector &)
{
  return CGAL_NTS is_zero(v.x()) && CGAL_NTS is_zero(v.y());
}

template < class R >
inline
bool
operator==(const Null_vector &n, const VectorC2<R> &v)
{
  return v == n;
}

template < class R >
inline
bool
operator!=(const VectorC2<R> &v, const Null_vector &n)
{
  return !(v == n);
}

template < class R >
inline
bool
operator!=(const Null_vector &n, const VectorC2<R> &v)
{
  return !(v == n);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_VECTOR_2_H
