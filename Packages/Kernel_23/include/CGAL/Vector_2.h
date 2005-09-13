// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// Author(s)     : Andreas Fabri, Stefan Schirra

#ifndef CGAL_VECTOR_2_H
#define CGAL_VECTOR_2_H

#include <CGAL/Origin.h>
#include <CGAL/Kernel/mpl.h>

CGAL_BEGIN_NAMESPACE

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
  typedef typename R_::Kernel_base::Vector_2  RVector_2;

public:

  typedef RVector_2 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_                        R;

  Vector_2() {}

  Vector_2(const Point_2& a, const Point_2& b)
      : RVector_2(typename R::Construct_vector_2()(a, b).rep()) {}

  Vector_2(const RVector_2& v) : RVector_2(v) {}

  Vector_2(const Segment_2 &s) : RVector_2(typename R::Construct_vector_2()(s).rep()) {}

  Vector_2(const Ray_2 &r) : RVector_2(typename R::Construct_vector_2()(r).rep()) {}

  Vector_2(const Line_2 &l) : RVector_2(typename R::Construct_vector_2()(l).rep()) {}

  Vector_2(const Null_vector &v) : RVector_2(typename R::Construct_vector_2()(v).rep()) {}

  Vector_2(const RT &x, const RT &y)
    : RVector_2(typename R::Construct_vector_2()(x,y).rep())
  {}

  Vector_2(const RT &x, const RT &y, const RT &w)
    : RVector_2(typename R::Construct_vector_2()(x,y,w).rep())
  {}


  typename Qualified_result_of<typename R::Compute_x_2,Vector_2>::type
  x() const
  {
    return R().compute_x_2_object()(*this);
  }

  
  typename Qualified_result_of<typename R::Compute_y_2,Vector_2>::type
  y() const
  {
    return R().compute_y_2_object()(*this);
  }

  
  typename Qualified_result_of<typename R::Compute_y_2,Vector_2>::type
  cartesian(int i) const
  {
    CGAL_kernel_precondition( (i == 0) || (i == 1) );
    return (i==0) ?  x() : y();
  }

  typename Qualified_result_of<typename R::Compute_x_2,Vector_2>::type
  operator[](int i) const
  {
      return cartesian(i);
  }

  typename Qualified_result_of<typename R::Compute_hx_2,Vector_2>::type
  hx() const
  {
    return R().compute_hx_2_object()(*this);
  }

  
  typename Qualified_result_of<typename R::Compute_hy_2,Vector_2>::type
  hy() const
  {
    return R().compute_hy_2_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hw_2,Vector_2>::type
  hw() const
  {
    return R().compute_hw_2_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hx_2,Vector_2>::type
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

  Vector_2 operator+(const Vector_2& v) const
  {
    return R().construct_sum_of_vectors_2_object()(*this,v);
  }

  Vector_2 operator/(const RT& c) const
  {
   return R().construct_divided_vector_2_object()(*this,c);
  }

  Vector_2 operator/(const typename First_if_different<FT,RT>::Type & c) const
  {
   return R().construct_divided_vector_2_object()(*this,c);
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



#ifndef CGAL_NO_OSTREAM_INSERT_VECTOR_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Vector_2<R> &v)
{
  return os << v.rep();
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTOR_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTOR_2
template < class R >
std::istream &
operator>>(std::istream &is, Vector_2<R>& v)
{
  return is >> v.rep();
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTOR_2

CGAL_END_NAMESPACE

#endif // CGAL_VECTOR_2_H
