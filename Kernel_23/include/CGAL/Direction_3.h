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

#ifndef CGAL_DIRECTION_3_H
#define CGAL_DIRECTION_3_H

#include <CGAL/assertions.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/IO/io.h>

namespace CGAL {

template <class R_>
class Direction_3 : public R_::Kernel_base::Direction_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Line_3                Line_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;

  typedef Direction_3                        Self;
  static_assert(std::is_same<Self, typename R_::Direction_3>::value);

public:

  typedef Dimension_tag<3>  Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  typedef typename R_::Kernel_base::Direction_3 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  Direction_3() {}

  Direction_3(const Rep& d)
    : Rep(d) {}

  Direction_3(Rep&& d)
    : Rep(std::move(d)) {}

  explicit Direction_3(const Vector_3& v)
    : Rep(typename R::Construct_direction_3()(Return_base_tag(), v)) {}

  explicit Direction_3(const Line_3& l)
    : Rep(typename R::Construct_direction_3()(Return_base_tag(), l)) {}

  explicit Direction_3(const Ray_3& r)
    : Rep(typename R::Construct_direction_3()(Return_base_tag(), r)) {}

  explicit Direction_3(const Segment_3& s)
    : Rep(typename R::Construct_direction_3()(Return_base_tag(), s)) {}

  Direction_3(const RT& hx, const RT& hy, const RT& hz)
    : Rep(typename R::Construct_direction_3()(Return_base_tag(), hx, hy, hz)) {}

  Direction_3 transform(const Aff_transformation_3 &t) const
  {
    return t.transform(*this);
  }

  Direction_3
  operator-() const
  {
    return R().construct_opposite_direction_3_object()(*this);
  }

  Vector_3 to_vector() const
  {
    return R().construct_vector_3_object()(*this);
  }

  Vector_3 vector() const { return to_vector(); }


  decltype(auto)
  dx() const
  {
    return R().compute_dx_3_object()(*this);
  }

  decltype(auto)
  dy() const
  {
    return R().compute_dy_3_object()(*this);
  }

  decltype(auto)
  dz() const
  {
    return R().compute_dz_3_object()(*this);
  }

  decltype(auto)
  delta(int i) const
  {
    CGAL_kernel_precondition( i >= 0 && i <= 2 );
    if (i==0) return dx();
    if (i==1) return dy();
    return dz();
  }

};


template <class R >
std::ostream&
insert(std::ostream& os, const Direction_3<R>& d, const Cartesian_tag&)
{
  typename R::Vector_3 v = d.to_vector();
  switch(IO::get_mode(os)) {
    case IO::ASCII :
      return os << v.x() << ' ' << v.y()  << ' ' << v.z();
    case IO::BINARY :
      write(os, v.x());
      write(os, v.y());
      write(os, v.z());
      return os;
    default:
      os << "DirectionC3(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
      return os;
  }
}

template <class R >
std::ostream&
insert(std::ostream& os, const Direction_3<R>& d, const Homogeneous_tag&)
{
  switch(IO::get_mode(os))
  {
    case IO::ASCII :
        return os << d.dx() << ' ' << d.dy() << ' ' << d.dz();
    case IO::BINARY :
        write(os, d.dx());
        write(os, d.dy());
        write(os, d.dz());
        return os;
    default:
        return os << "DirectionH3(" << d.dx() << ", "
                                    << d.dy() << ", "
                                    << d.dz() << ')';
  }
}

template < class R >
std::ostream&
operator<<(std::ostream& os, const Direction_3<R>& d)
{
  return insert(os, d, typename R::Kernel_tag() );
}


template <class R >
std::istream&
extract(std::istream& is, Direction_3<R>& d, const Cartesian_tag&)
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
      d = Direction_3<R>(x, y, z);
  return is;
}

template <class R >
std::istream&
extract(std::istream& is, Direction_3<R>& d, const Homogeneous_tag&)
{
  typename R::RT x, y, z;
  switch(IO::get_mode(is))
  {
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
    d = Direction_3<R>(x, y, z);
  return is;
}

template < class R >
std::istream&
operator>>(std::istream& is, Direction_3<R>& d)
{
  return extract(is, d, typename R::Kernel_tag() );
}

} //namespace CGAL

#endif // CGAL_DIRECTION_3_H
