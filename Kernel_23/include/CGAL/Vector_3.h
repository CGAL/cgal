// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
#include <CGAL/Kernel/mpl.h>
#include <CGAL/representation_tags.h>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R_>
class Vector_3 : public R_::Kernel_base::Vector_3
{
  typedef typename R_::RT                    RT;
// http://doc.cgal.org/latest/Manual/devman_code_format.html#secprogramming_conventions
  typedef typename R_::FT                    FT_;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Line_3                Line_3;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;

  typedef Vector_3                            Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Vector_3>::value));

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
  Vector_3(const T1 &x, const T2 &y, const T3 &z)
    : Rep(typename R::Construct_vector_3()(Return_base_tag(), x, y, z)) {}

  Vector_3(const RT& x, const RT& y, const RT& z, const RT& w)
    : Rep(typename R::Construct_vector_3()(Return_base_tag(), x, y, z, w)) {}

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

  Vector_3 operator/(const typename First_if_different<FT_,RT>::Type & c) const
  {
   return R().construct_divided_vector_3_object()(*this,c);
  }

  typename cpp11::result_of<typename R::Compute_x_3(Vector_3)>::type
  x() const
  {
    return R().compute_x_3_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_y_3(Vector_3)>::type
  y() const
  {
    return R().compute_y_3_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_z_3(Vector_3)>::type
  z() const
  {
    return R().compute_z_3_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_hx_3(Vector_3)>::type
  hx() const
  {
    return R().compute_hx_3_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_hy_3(Vector_3)>::type
  hy() const
  {
    return R().compute_hy_3_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_hz_3(Vector_3)>::type
  hz() const
  {
    return R().compute_hz_3_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_hw_3(Vector_3)>::type
  hw() const
  {
    return R().compute_hw_3_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_x_3(Vector_3)>::type
  cartesian(int i) const
  {
    CGAL_kernel_precondition( (i == 0) || (i == 1) || (i == 2) );
    if (i==0) return x();
    if (i==1) return y();
    return z();
  }

  typename cpp11::result_of<typename R::Compute_hw_3(Vector_3)>::type
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

  typename cpp11::result_of<typename R::Compute_x_3(Vector_3)>::type
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

  typename cpp11::result_of<typename R::Compute_squared_length_3(Vector_3)>::type
  squared_length() const
  {
    return R().compute_squared_length_3_object()(*this);
  }

};


template <class R >
std::ostream&
insert(std::ostream& os, const Vector_3<R>& v, const Cartesian_tag&) 
{
  switch(os.iword(IO::mode)) {
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
  switch(os.iword(IO::mode))
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
  typename R::FT x, y, z;
  switch(is.iword(IO::mode)) {
    case IO::ASCII :
      is >> iformat(x) >> iformat(y) >> iformat(z);
      break;
    case IO::BINARY :
      read(is, x);
      read(is, y);
      read(is, z);
      break;
    default:
      std::cerr << "" << std::endl;
      std::cerr << "Stream must be in ascii or binary mode" << std::endl;
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
  switch(is.iword(IO::mode))
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
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
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
