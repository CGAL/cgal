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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri, Stefan Schirra

#ifndef CGAL_VECTOR_3_H
#define CGAL_VECTOR_3_H

#include <CGAL/Origin.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Vector_3 : public R_::Kernel_base::Vector_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::FT                    FT;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Line_3                Line_3;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;
public:

  typedef typename R_::Kernel_base::Vector_3         Rep;

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

  Vector_3(const Point_3& a, const Point_3& b)
    : Rep(typename R::Construct_vector_3()(a, b).rep()) {}

  Vector_3(const Segment_3& s)
    : Rep(typename R::Construct_vector_3()(s).rep()) {}

  Vector_3(const Ray_3& r)
    : Rep(typename R::Construct_vector_3()(r).rep()) {}

  Vector_3(const Line_3& l)
    : Rep(typename R::Construct_vector_3()(l).rep()) {}

  Vector_3(const Null_vector& v)
    : Rep(typename R::Construct_vector_3()(v).rep()) {}

  Vector_3(const RT& x, const RT& y, const RT& z)
    : Rep(typename R::Construct_vector_3()(x, y, z).rep()) {}

  Vector_3(const RT& x, const RT& y, const RT& z, const RT& w)
    : Rep(typename R::Construct_vector_3()(x, y, z, w).rep()) {}

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

  Vector_3 operator+(const Vector_3& v) const
  {
    return R().construct_sum_of_vectors_3_object()(*this,v);
  }

  Vector_3 operator/(const RT& c) const
  {
   return R().construct_divided_vector_3_object()(*this,c);
  }

  Vector_3 operator/(const typename First_if_different<FT,RT>::Type & c) const
  {
   return R().construct_divided_vector_3_object()(*this,c);
  }

  typename Qualified_result_of<typename R::Compute_x_3, Vector_3>::type
  x() const
  {
    return R().compute_x_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_y_3, Vector_3>::type
  y() const
  {
    return R().compute_y_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_z_3, Vector_3>::type
  z() const
  {
    return R().compute_z_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hx_3, Vector_3>::type
  hx() const
  {
    return R().compute_hx_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hy_3, Vector_3>::type
  hy() const
  {
    return R().compute_hy_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hz_3, Vector_3>::type
  hz() const
  {
    return R().compute_hz_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hw_3, Vector_3>::type
  hw() const
  {
    return R().compute_hw_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_x_3, Vector_3>::type
  cartesian(int i) const
  {
    CGAL_kernel_precondition( (i == 0) || (i == 1) || (i == 2) );
    if (i==0) return x();
    if (i==1) return y();
    return z();
  }

  typename Qualified_result_of<typename R::Compute_hw_3, Vector_3>::type
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

  typename Qualified_result_of<typename R::Compute_x_3, Vector_3>::type
  operator[](int i) const
  {
      return cartesian(i);
  }

  typename Qualified_result_of<typename R::Compute_squared_length_3, Vector_3>::type
  squared_length() const
  {
    return R().compute_squared_length_3_object()(*this);
  }

};

#ifndef CGAL_NO_OSTREAM_INSERT_VECTOR_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Vector_3<R>& v)
{
  typedef typename  R::Kernel_base::Vector_3  Rep;
  return os << static_cast<const Rep&>(v);
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTOR_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTOR_3
template < class R >
std::istream&
operator>>(std::istream& is, Vector_3<R>& p)
{
  typedef typename  R::Kernel_base::Vector_3  Rep;
  return is >> static_cast<Rep&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTOR_3

CGAL_END_NAMESPACE

#endif // CGAL_VECTOR_3_H
