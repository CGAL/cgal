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

#ifndef CGAL_VECTOR_2_H
#define CGAL_VECTOR_2_H

#include <CGAL/Origin.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/result_of.h>
#include <CGAL/IO/io.h>

namespace CGAL {

template <class R_>
class Vector_2 : public R_::Kernel_base::Vector_2
{
  typedef typename R_::RT             RT;
  typedef typename R_::FT             FT;
  typedef typename R_::Segment_2      Segment_2;
  typedef typename R_::Ray_2          Ray_2;
  typedef typename R_::Line_2         Line_2;
  typedef typename R_::Point_2        Point_2;
  typedef typename R_::Direction_2    Direction_2;
  typedef typename R_::Aff_transformation_2  Aff_transformation_2;
  typedef typename R_::Kernel_base::Vector_2  RVector_2;

  typedef Vector_2                    Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Vector_2>::value));

public:

  typedef Dimension_tag<2>  Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  typedef RVector_2 Rep;
  typedef typename R_::Cartesian_const_iterator_2 Cartesian_const_iterator;

  const Rep& rep() const noexcept
  {
    return *this;
  }

  Rep& rep() noexcept
  {
    return *this;
  }

  typedef  R_                        R;

  Vector_2() {}

  Vector_2(const RVector_2& v)
      : RVector_2(v) {}

  Vector_2(const Point_2& a, const Point_2& b)
      : RVector_2(typename R::Construct_vector_2()(Return_base_tag(), a, b)) {}

  explicit Vector_2(const Segment_2 &s)
      : RVector_2(typename R::Construct_vector_2()(Return_base_tag(), s)) {}

  explicit Vector_2(const Ray_2 &r)
      : RVector_2(typename R::Construct_vector_2()(Return_base_tag(), r)) {}

  explicit Vector_2(const Line_2 &l)
      : RVector_2(typename R::Construct_vector_2()(Return_base_tag(), l)) {}

  Vector_2(const Null_vector &v)
      : RVector_2(typename R::Construct_vector_2()(Return_base_tag(), v)) {}

  template < typename T1, typename T2 >
  Vector_2(const T1 &x, const T2 &y)
      : RVector_2(typename R::Construct_vector_2()(Return_base_tag(), x,y)) {}

  Vector_2(const RT &x, const RT &y, const RT &w)
      : RVector_2(typename R::Construct_vector_2()(Return_base_tag(), x,y,w)) {}

  friend void swap(Self& a, Self& b)
#ifdef __cpp_lib_is_swappable
    noexcept(std::is_nothrow_swappable_v<Rep>)
#endif
  {
    using std::swap;
    swap(a.rep(), b.rep());
  }


  typename cpp11::result_of<typename R::Compute_x_2(Vector_2)>::type
  x() const
  {
    return R().compute_x_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_y_2(Vector_2)>::type
  y() const
  {
    return R().compute_y_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_y_2(Vector_2)>::type
  cartesian(int i) const
  {
    CGAL_kernel_precondition( (i == 0) || (i == 1) );
    return (i==0) ?  x() : y();
  }

  typename cpp11::result_of<typename R::Compute_x_2(Vector_2)>::type
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

  typename cpp11::result_of<typename R::Compute_hx_2(Vector_2)>::type
  hx() const
  {
    return R().compute_hx_2_object()(*this);
  }


  typename cpp11::result_of<typename R::Compute_hy_2(Vector_2)>::type
  hy() const
  {
    return R().compute_hy_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_hw_2(Vector_2)>::type
  hw() const
  {
    return R().compute_hw_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_hx_2(Vector_2)>::type
  homogeneous(int i) const
  {
    CGAL_kernel_precondition( (i >= 0) || (i <= 2) );
    return (i==0) ?  hx() : (i==1)? hy() : hw();
  }

  int dimension() const
  {
      return 2;
  }

  Vector_2 operator-() const
  {
    return R().construct_opposite_vector_2_object()(*this);
  }

  Vector_2 operator-(const Vector_2& v) const
  {
    return R().construct_difference_of_vectors_2_object()(*this,v);
  }

  Vector_2& operator-=(const Vector_2& v)
  {
    *this = R().construct_difference_of_vectors_2_object()(*this,v);
    return *this;
  }

  Vector_2 operator+(const Vector_2& v) const
  {
    return R().construct_sum_of_vectors_2_object()(*this,v);
  }

  Vector_2& operator+=(const Vector_2& v)
  {
    *this = R().construct_sum_of_vectors_2_object()(*this,v);
    return *this;
  }

  Vector_2 operator/(const RT& c) const
  {
   return R().construct_divided_vector_2_object()(*this,c);
  }

  Vector_2& operator/=(const RT& c)
  {
    *this = R().construct_divided_vector_2_object()(*this,c);
    return *this;
  }

  Vector_2 operator/(const typename First_if_different<FT,RT>::Type & c) const
  {
   return R().construct_divided_vector_2_object()(*this,c);
  }

  Vector_2& operator/=(const typename First_if_different<FT,RT>::Type & c)
  {
    *this = R().construct_divided_vector_2_object()(*this,c);
    return *this;
  }

  Vector_2& operator*=(const RT& c)
  {
    *this = R().construct_scaled_vector_2_object()(*this,c);
    return *this;
  }

  Vector_2& operator*=(const typename First_if_different<FT,RT>::Type & c)
  {
    *this = R().construct_scaled_vector_2_object()(*this,c);
    return *this;
  }

  FT squared_length() const
  {
    return R().compute_squared_length_2_object()(*this);
  }


  Direction_2 direction() const
  {
    return R().construct_direction_2_object()(*this);
  }

  Vector_2 perpendicular(const Orientation &o) const
  {
    return R().construct_perpendicular_vector_2_object()(*this,o);
  }

  Vector_2 transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }

};


template < class R >
inline
bool
operator==(const Vector_2<R> &v, const Null_vector &n)
{
  return R().equal_2_object()(v, n);
}

template < class R >
inline
bool
operator==(const Null_vector &n, const Vector_2<R> &v)
{
  return v == n;
}

template < class R >
inline
bool
operator!=(const Vector_2<R> &v, const Null_vector &n)
{
  return !(v == n);
}

template < class R >
inline
bool
operator!=(const Null_vector &n, const Vector_2<R> &v)
{
  return !(v == n);
}


template <class R >
std::ostream&
insert(std::ostream& os, const Vector_2<R>& v, const Cartesian_tag&)
{
    switch(get_mode(os)) {
    case IO::ASCII :
        return os << v.x() << ' ' << v.y();
    case IO::BINARY :
        write(os, v.x());
        write(os, v.y());
        return os;
    default:
        return os << "VectorC2(" << v.x() << ", " << v.y() << ')';
    }
}

template <class R >
std::ostream&
insert(std::ostream& os, const Vector_2<R>& v, const Homogeneous_tag&)
{
  switch(get_mode(os))
  {
    case IO::ASCII :
        return os << v.hx() << ' ' << v.hy() << ' ' << v.hw();
    case IO::BINARY :
        write(os, v.hx());
        write(os, v.hy());
        write(os, v.hw());
        return os;
    default:
        return os << "VectorH2(" << v.hx() << ", "
                                 << v.hy() << ", "
                                 << v.hw() << ')';
  }
}

template < class R >
std::ostream&
operator<<(std::ostream& os, const Vector_2<R>& v)
{
  return insert(os, v, typename R::Kernel_tag() );
}



template <class R >
std::istream&
extract(std::istream& is, Vector_2<R>& v, const Cartesian_tag&)
{
  typename R::FT x(0), y(0);
    switch(get_mode(is)) {
    case IO::ASCII :
        is >> iformat(x) >> iformat(y);
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        is.setstate(std::ios::failbit);
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
        v = Vector_2<R>(x, y);
    return is;
}


template <class R >
std::istream&
extract(std::istream& is, Vector_2<R>& v, const Homogeneous_tag&)
{
  typename R::RT hx, hy, hw;
  switch(get_mode(is))
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
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  v = Vector_2<R>(hx, hy, hw);
  return is;
}

template < class R >
std::istream&
operator>>(std::istream& is, Vector_2<R>& v)
{
  return extract(is, v, typename R::Kernel_tag() );
}

} //namespace CGAL

#endif // CGAL_VECTOR_2_H
