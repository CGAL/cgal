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
// Author(s)     : Andreas Fabri and Herve Bronnimann

#ifndef CGAL_CARTESIAN_POINT_3_H
#define CGAL_CARTESIAN_POINT_3_H

#include <CGAL/Origin.h>

namespace CGAL {

template < class R_ >
class PointC3
{
  typedef PointC3<R_>                       Self;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  // We do not use reference counting here as it is done at the Vector_3 level.
  Vector_3 base;

public:
  typedef typename Vector_3::Cartesian_const_iterator Cartesian_const_iterator;
  typedef R_                                R;
  typedef typename R_::FT                   FT;

  PointC3() {}

  PointC3(const Origin &)
    : base(NULL_VECTOR) {}

  PointC3(const FT &x, const FT &y, const FT &z)
    : base(x, y, z) {}

  PointC3(const FT &x, const FT &y, const FT &z, const FT &w)
    : base(x, y, z, w) {}

  friend void swap(Self& a, Self& b)
#ifdef __cpp_lib_is_swappable
    noexcept(std::is_nothrow_swappable_v<Vector_3>)
#endif
  {
    using std::swap;
    swap(a.base, b.base);
  }

  const FT & x() const
  {
      return base.x();
  }
  const FT & y() const
  {
      return base.y();
  }
  const FT & z() const
  {
      return base.z();
  }

  const FT & hx() const
  {
      return base.hx();
  }
  const FT & hy() const
  {
      return base.hy();
  }
  const FT & hz() const
  {
      return base.hz();
  }
  const FT & hw() const
  {
      return base.hw();
  }

  const FT & cartesian(int i) const;
  const FT & operator[](int i) const;
  const FT & homogeneous(int i) const;

  Cartesian_const_iterator cartesian_begin() const
  {
    return base.cartesian_begin();
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return base.cartesian_end();
  }

  int dimension() const
  {
      return base.dimension();
  }

  Point_3 transform(const Aff_transformation_3 &t) const
  {
    return t.transform(*this);
  }
};

template < class R >
inline
const typename PointC3<R>::FT &
PointC3<R>::cartesian(int i) const
{
  return base.cartesian(i);
}

template < class R >
inline
const typename PointC3<R>::FT &
PointC3<R>::operator[](int i) const
{
  return base[i];
}

template < class R >
inline
const typename PointC3<R>::FT &
PointC3<R>::homogeneous(int i) const
{
  return base.homogeneous(i);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_POINT_3_H
