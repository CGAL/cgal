// Copyright (c) 1999
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
// Author(s)     : Andreas Fabri, Stefan Schirra

#ifndef CGAL_POINT_2_H
#define CGAL_POINT_2_H

#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/assertions.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Dimension.h>

#include <type_traits>

namespace CGAL {

template <class R_>
class Point_2 : public R_::Kernel_base::Point_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::FT                    FT;
  typedef typename R_::Aff_transformation_2  Aff_transformation_2;
  typedef typename R_::Kernel_base::Point_2  RPoint_2;

  typedef Point_2                            Self;
  static_assert(std::is_same<Self, typename R_::Point_2>::value);

public:

  typedef Dimension_tag<2>  Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  typedef typename R_::Weighted_point_2 Weighted_point_2;
  typedef RPoint_2 Rep;
  typedef typename R_::Cartesian_const_iterator_2 Cartesian_const_iterator;

  const Rep& rep() const noexcept
  {
    return *this;
  }

  Rep& rep() noexcept
  {
    return *this;
  }

  typedef  R_   R;

  Point_2() {}

  Point_2(const Origin& o)
    : RPoint_2(typename R::Construct_point_2()(Return_base_tag(), o))
  {}

  Point_2(const RPoint_2& p)
    : RPoint_2(p)
  {}

  Point_2(RPoint_2&& p)
    : RPoint_2(std::move(p))
  {}

  explicit
  Point_2(const Weighted_point_2& wp)
    : Rep(wp.point())
  {}

  template < typename T1, typename T2>
  Point_2(T1&& x, T2&& y)
    : Rep(typename R::Construct_point_2()(Return_base_tag(),
                                          std::forward<T1>(x),
                                          std::forward<T2>(y)))
  {}

  Point_2(const RT& hx, const RT& hy, const RT& hw)
    : RPoint_2(typename R::Construct_point_2()(Return_base_tag(), hx, hy, hw))
  {}

  friend void swap(Self& a, Self& b)
#if !defined(__INTEL_COMPILER) && defined(__cpp_lib_is_swappable)
    noexcept(std::is_nothrow_swappable_v<Rep>)
#endif
  {
    using std::swap;
    swap(a.rep(), b.rep());
  }

  decltype(auto)
  x() const
  {
    return typename R::Compute_x_2()(*this);
  }

  decltype(auto)
  y() const
  {
    return typename R::Compute_y_2()(*this);
  }

  decltype(auto)
  cartesian(int i) const
  {
    CGAL_kernel_precondition( (i == 0) || (i == 1) );
    return (i==0) ?  x() : y();
  }

  decltype(auto)
  operator[](int i) const
  {
      return cartesian(i);
  }

  Cartesian_const_iterator cartesian_begin() const
  {
    return typename R::Construct_cartesian_const_iterator_2()(*this);
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return typename R::Construct_cartesian_const_iterator_2()(*this,2);
  }


  decltype(auto)
  hx() const
  {
    return typename R::Compute_hx_2()(*this);
  }

  decltype(auto)
  hy() const
  {
    return typename R::Compute_hy_2()(*this);
  }

  decltype(auto)
  hw() const
  {
    return typename R::Compute_hw_2()(*this);
  }

  int dimension() const
  {
      return 2;
  }

  decltype(auto)
  homogeneous(int i) const
  {
    CGAL_kernel_precondition( (i >= 0) || (i <= 2) );
    return (i==0) ?  hx() : (i==1)? hy() : hw();
  }

  Bbox_2 bbox() const
  {
    return R().construct_bbox_2_object()(*this);
  }

  Point_2 transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }

  Point_2<R_>&
  operator+=(const typename R::Vector_2 &v)
  {
    *this = R().construct_translated_point_2_object()(*this, v);
    return *this;
  }

  Point_2<R_>&
  operator-=(const typename R::Vector_2 &v)
  {
    *this = R().construct_translated_point_2_object()(*this,
                  R().construct_opposite_vector_2_object()(v));
    return *this;
  }

};


template <class R >
std::ostream&
insert(std::ostream& os, const Point_2<R>& p,const Cartesian_tag&)
{
    switch(IO::get_mode(os)) {
    case IO::ASCII :
        return os << p.x() << ' ' << p.y();
    case IO::BINARY :
        write(os, p.x());
        write(os, p.y());
        return os;
    default:
        return os << "PointC2(" << p.x() << ", " << p.y() << ')';
    }
}

template <class R >
std::ostream&
insert(std::ostream& os, const Point_2<R>& p,const Homogeneous_tag&)
{
  switch(IO::get_mode(os))
  {
    case IO::ASCII :
        return os << p.hx() << ' ' << p.hy() << ' ' << p.hw();
    case IO::BINARY :
        write(os, p.hx());
        write(os, p.hy());
        write(os, p.hw());
        return os;
    default:
        return os << "PointH2(" << p.hx() << ", "
                                << p.hy() << ", "
                                << p.hw() << ')';
  }
}

template < class R >
std::ostream&
operator<<(std::ostream& os, const Point_2<R>& p)
{
  return insert(os, p, typename R::Kernel_tag() );
}


template <class R >
std::istream&
extract(std::istream& is, Point_2<R>& p, const Cartesian_tag&)
{
  typename R::FT x(0), y(0);
    switch(IO::get_mode(is)) {
    case IO::ASCII :
        is >> IO::iformat(x) >> IO::iformat(y);
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        is.setstate(std::ios::failbit);
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ASCII or binary mode" << std::endl;
        break;
    }
    if (is)
      p = Point_2<R>(std::move(x), std::move(y));
    return is;
}


template <class R >
std::istream&
extract(std::istream& is, Point_2<R>& p, const Homogeneous_tag&)
{
  typename R::RT hx, hy, hw;
  switch(IO::get_mode(is))
  {
    case IO::ASCII :
        is >> hx >> hy >> hw;
        break;
    case IO::BINARY :
        read(is, hx);
        read(is, hy);
        read(is, hw);
        break;
    default:
        is.setstate(std::ios::failbit);
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ASCII or binary mode" << std::endl;
        break;
  }
  if (is)
    p = Point_2<R>(hx, hy, hw);
  return is;
}

template < class R >
std::istream&
operator>>(std::istream& is, Point_2<R>& p)
{
  return extract(is, p, typename R::Kernel_tag() );
}

} //namespace CGAL

#endif // CGAL_POINT_2_H
