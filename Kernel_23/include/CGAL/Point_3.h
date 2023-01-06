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

#ifndef CGAL_POINT_3_H
#define CGAL_POINT_3_H

#include <CGAL/Origin.h>
#include <CGAL/representation_tags.h>
#include <CGAL/assertions.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R_>
class Point_3 : public R_::Kernel_base::Point_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;

  typedef Point_3                            Self;
  CGAL_static_assertion((std::is_same<Self, typename R_::Point_3>::value));

public:

  typedef Dimension_tag<3>  Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  typedef typename R_::Weighted_point_3 Weighted_point_3;
  typedef typename R_::Kernel_base::Point_3  Rep;
  typedef typename R_::Cartesian_const_iterator_3 Cartesian_const_iterator;

  const Rep& rep() const noexcept
  {
    return *this;
  }

  Rep& rep() noexcept
  {
    return *this;
  }

  typedef          R_                       R;

  Point_3() {}

  Point_3(const Origin& o)
    : Rep(typename R::Construct_point_3()(Return_base_tag(), o))
  {}

  Point_3(const Rep& p)
      : Rep(p) {}

  Point_3(Rep&& p)
      : Rep(std::move(p)) {}

  explicit
  Point_3(const Weighted_point_3& wp)
    : Rep(wp.point())
  {}

  template < typename T1, typename T2, typename T3 >
  Point_3(const T1& x, const T2& y, const T3& z)
    : Rep(typename R::Construct_point_3()(Return_base_tag(), x, y, z))
  {}

  Point_3(const RT& hx, const RT& hy, const RT& hz, const RT& hw)
    : Rep(typename R::Construct_point_3()(Return_base_tag(), hx, hy, hz, hw))
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
    return typename R::Compute_x_3()(*this);
  }

  decltype(auto)
  y() const
  {
    return typename R::Compute_y_3()(*this);
  }

  decltype(auto)
  z() const
  {
    return typename R::Compute_z_3()(*this);
  }

  decltype(auto)
  hx() const
  {
    return R().compute_hx_3_object()(*this);
  }

  decltype(auto)
  hy() const
  {
    return R().compute_hy_3_object()(*this);
  }

  decltype(auto)
  hz() const
  {
    return R().compute_hz_3_object()(*this);
  }

  decltype(auto)
  hw() const
  {
    return R().compute_hw_3_object()(*this);
  }

  decltype(auto)
  cartesian(int i) const
  {
    CGAL_kernel_precondition( (i == 0) || (i == 1) || (i == 2) );
    if (i==0) return x();
    if (i==1) return y();
    return z();
  }

  RT
  homogeneous(int i) const
  {
    CGAL_kernel_precondition( (i >= 0) || (i <= 3) );
    if (i==0) return hx();
    if (i==1) return hy();
    if (i==2) return hz();
    return hw();
  }

  decltype(auto)
  operator[](int i) const
  {
      return cartesian(i);
  }

  Cartesian_const_iterator cartesian_begin() const
  {
    return typename R::Construct_cartesian_const_iterator_3()(*this);
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return typename R::Construct_cartesian_const_iterator_3()(*this,3);
  }

  int dimension() const
  {
      return 3;
  }

  Bbox_3 bbox() const
  {
    return R().construct_bbox_3_object()(*this);
  }

  Point_3 transform(const Aff_transformation_3 &t) const
  {
    return t.transform(*this);
  }

  Point_3<R_>&
  operator+=(const Vector_3 &v)
  {
    *this = R().construct_translated_point_3_object()(*this, v);
    return *this;
  }

  Point_3<R_>&
  operator-=(const Vector_3 &v)
  {
    *this = R().construct_translated_point_3_object()(*this,
                  R().construct_opposite_vector_3_object()(v));
    return *this;
  }

};

template <class R>
inline
bool
operator==(const Origin& o, const Point_3<R>& p)
{ return p == o; }

template <class R>
inline
bool
operator!=(const Origin& o, const Point_3<R>& p)
{ return p != o; }


template <class R >
std::ostream&
insert(std::ostream& os, const Point_3<R>& p,const Cartesian_tag&)
{
    switch(IO::get_mode(os)) {
    case IO::ASCII :
        return os << p.x() << ' ' << p.y() << ' ' << p.z();
    case IO::BINARY :
        write(os, p.x());
        write(os, p.y());
        write(os, p.z());
        return os;
    default:
        return os << "PointC3(" << p.x() << ", " << p.y()
                                         << ", " << p.z() << ')';
    }
}

template <class R >
std::ostream&
insert(std::ostream& os, const Point_3<R>& p,const Homogeneous_tag&)
{
  switch(IO::get_mode(os))
  {
    case IO::ASCII :
        return os << p.hx() << ' ' << p.hy() << ' ' << p.hz() << ' ' << p.hw();
    case IO::BINARY :
        write(os, p.hx());
        write(os, p.hy());
        write(os, p.hz());
        write(os, p.hw());
        return os;
    default:
        return os << "PointH3(" << p.hx() << ", "
                                << p.hy() << ", "
                                << p.hz() << ", "
                                << p.hw() << ')';
  }
}

template < class R >
std::ostream&
operator<<(std::ostream& os, const Point_3<R>& p)
{
  return insert(os, p, typename R::Kernel_tag() );
}


template <class R >
std::istream&
extract(std::istream& is, Point_3<R>& p, const Cartesian_tag&)
{
  typename R::FT x(0), y(0), z(0);
    switch(IO::get_mode(is)) {
    case IO::ASCII :
        is >> IO::iformat(x) >> IO::iformat(y) >> IO::iformat(z);
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        read(is, z);
        break;
    default:
        is.setstate(std::ios::failbit);
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ASCII or binary mode" << std::endl;
        break;
    }
    if (is)
      p = Point_3<R>(x, y, z);
    return is;
}


template <class R >
std::istream&
extract(std::istream& is, Point_3<R>& p, const Homogeneous_tag&)
{
  typename R::RT hx, hy, hz, hw;
  switch(IO::get_mode(is))
  {
    case IO::ASCII :
        is >> hx >> hy >> hz >> hw;
        break;
    case IO::BINARY :
        read(is, hx);
        read(is, hy);
        read(is, hz);
        read(is, hw);
        break;
    default:
        is.setstate(std::ios::failbit);
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ASCII or binary mode" << std::endl;
        break;
  }
  if (is)
    p = Point_3<R>(hx, hy, hz, hw);
  return is;
}

template < class R >
std::istream&
operator>>(std::istream& is, Point_3<R>& p)
{
  return extract(is, p, typename R::Kernel_tag() );
}

} //namespace CGAL

#endif // CGAL_POINT_3_H
