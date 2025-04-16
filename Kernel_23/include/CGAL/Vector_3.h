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

#ifndef CGAL_VECTOR_3_H
#define CGAL_VECTOR_3_H

#include <CGAL/Origin.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/representation_tags.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/assertions.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>
#include <CGAL/IO/io.h>

namespace CGAL {

template <class R_>
class Vector_3 : public R_::Kernel_base::Vector_3
{
  typedef typename R_::RT                    RT;
// https://doc.cgal.org/latest/Manual/devman_code_format.html#secprogramming_conventions
  typedef typename R_::FT                    FT_;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Line_3                Line_3;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;

  typedef Vector_3                            Self;
  static_assert(std::is_same<Self, typename R_::Vector_3>::value);

public:

  typedef Dimension_tag<3>  Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  typedef typename R_::Cartesian_const_iterator_3 Cartesian_const_iterator;
  typedef typename R_::Kernel_base::Vector_3      Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  Vector_3() {}

  Vector_3(const Rep& v)
      : Rep(v) {}

  Vector_3(Rep&& v)
      : Rep(std::move(v)) {}

  Vector_3(const Point_3& a, const Point_3& b)
    : Rep(typename R::Construct_vector_3()(Return_base_tag(), a, b)) {}

  explicit Vector_3(const Segment_3& s)
    : Rep(typename R::Construct_vector_3()(Return_base_tag(), s)) {}

  explicit Vector_3(const Ray_3& r)
    : Rep(typename R::Construct_vector_3()(Return_base_tag(), r)) {}

  explicit Vector_3(const Line_3& l)
    : Rep(typename R::Construct_vector_3()(Return_base_tag(), l)) {}

  Vector_3(const Null_vector& v)
    : Rep(typename R::Construct_vector_3()(Return_base_tag(), v)) {}

  template < typename T1, typename T2, typename T3 >
  Vector_3(T1&& x, T2&& y, T3&& z)
    : Rep(typename R::Construct_vector_3()(Return_base_tag(), std::forward<T1>(x),
                                                              std::forward<T2>(y),
                                                              std::forward<T3>(z)))
  {}

  Vector_3(const RT& x, const RT& y, const RT& z, const RT& w)
    : Rep(typename R::Construct_vector_3()(Return_base_tag(), x, y, z, w)) {}

  friend void swap(Self& a, Self& b)
#if !defined(__INTEL_COMPILER) && defined(__cpp_lib_is_swappable)
    noexcept(std::is_nothrow_swappable_v<Rep>)
#endif
  {
    using std::swap;
    swap(a.rep(), b.rep());
  }

  Direction_3 direction() const
  {
    return R().construct_direction_3_object()(*this);
  }

  Vector_3 transform(const Aff_transformation_3 &t) const
  {
    return t.transform(*this);
  }

  Vector_3 operator-() const
  {
    return R().construct_opposite_vector_3_object()(*this);
  }

  Vector_3 operator-(const Vector_3& v) const
  {
    return R().construct_difference_of_vectors_3_object()(*this,v);
  }

  Vector_3& operator-=(const Vector_3& v)
  {
    *this = R().construct_difference_of_vectors_3_object()(*this,v);
    return *this;
  }

  Vector_3 operator+(const Vector_3& v) const
  {
    return R().construct_sum_of_vectors_3_object()(*this,v);
  }

  Vector_3& operator+=(const Vector_3& v)
  {
    *this = R().construct_sum_of_vectors_3_object()(*this,v);
    return *this;
  }

  Vector_3 operator/(const RT& c) const
  {
   return R().construct_divided_vector_3_object()(*this,c);
  }

  Vector_3& operator/=(const RT& c)
  {
    *this = R().construct_divided_vector_3_object()(*this,c);
    return *this;
  }

  Vector_3 operator/(const typename First_if_different<FT_,RT>::Type & c) const
  {
   return R().construct_divided_vector_3_object()(*this,c);
  }

  Vector_3& operator/=(const typename First_if_different<FT_,RT>::Type & c)
  {
    *this = R().construct_divided_vector_3_object()(*this,c);
    return *this;
  }

  Vector_3& operator*=(const RT& c)
  {
    *this = R().construct_scaled_vector_3_object()(*this,c);
    return *this;
  }

  Vector_3& operator*=(const typename First_if_different<FT_,RT>::Type & c)
  {
    *this = R().construct_scaled_vector_3_object()(*this,c);
    return *this;
  }

  decltype(auto)
  x() const
  {
    return R().compute_x_3_object()(*this);
  }

  decltype(auto)
  y() const
  {
    return R().compute_y_3_object()(*this);
  }

  decltype(auto)
  z() const
  {
    return R().compute_z_3_object()(*this);
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

  decltype(auto)
  homogeneous(int i) const
  {
    CGAL_kernel_precondition( (i >= 0) || (i <= 3) );
    if (i==0) return hx();
    if (i==1) return hy();
    if (i==2) return hz();
    return hw();
  }

  int dimension() const
  {
      return 3;
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

  decltype(auto)
  squared_length() const
  {
    return R().compute_squared_length_3_object()(*this);
  }

};


template <class R >
std::ostream&
insert(std::ostream& os, const Vector_3<R>& v, const Cartesian_tag&)
{
  switch(IO::get_mode(os)) {
    case IO::ASCII :
      return os << v.x() << ' ' << v.y()  << ' ' << v.z();
    case IO::BINARY :
      write(os, v.x());
      write(os, v.y());
      write(os, v.z());
      return os;
    default:
      os << "VectorC3(" << v.x() << ", " << v.y() <<  ", " << v.z() << ")";
      return os;
  }
}

template <class R >
std::ostream&
insert(std::ostream& os, const Vector_3<R>& v, const Homogeneous_tag&)
{
  switch(IO::get_mode(os))
  {
    case IO::ASCII :
        return os << v.hx() << ' ' << v.hy() << ' ' << v.hz() << ' ' << v.hw();
    case IO::BINARY :
        write(os, v.hx());
        write(os, v.hy());
        write(os, v.hz());
        write(os, v.hw());
        return os;
    default:
        return os << "VectorH3(" << v.hx() << ", "
                                 << v.hy() << ", "
                                 << v.hz() << ", "
                                 << v.hw() << ')';
  }
}

template < class R >
std::ostream&
operator<<(std::ostream& os, const Vector_3<R>& v)
{
  return insert(os, v, typename R::Kernel_tag() );
}


template <class R >
std::istream&
extract(std::istream& is, Vector_3<R>& v, const Cartesian_tag&)
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
      v = Vector_3<R>(x, y, z);
  return is;
}

template <class R >
std::istream&
extract(std::istream& is, Vector_3<R>& v, const Homogeneous_tag&)
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
    v = Vector_3<R>(hx, hy, hz, hw);
  return is;
}

template < class R >
std::istream&
operator>>(std::istream& is, Vector_3<R>& v)
{
  return extract(is, v, typename R::Kernel_tag() );
}

} //namespace CGAL

#endif // CGAL_VECTOR_3_H
