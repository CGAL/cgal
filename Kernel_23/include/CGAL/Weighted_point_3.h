// Copyright (c) 2016, 2017
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
// Author(s)     : Mariette Yvinec, Sylvain Pion

#ifndef CGAL_WEIGHTED_POINT_3_H
#define CGAL_WEIGHTED_POINT_3_H

#include <CGAL/Origin.h>
#include <CGAL/representation_tags.h>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Dimension.h>
#include <CGAL/Point_3.h>

namespace CGAL {

template <class R_>
class Weighted_point_3 : public R_::Kernel_base::Weighted_point_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::FT                    FT;

  typedef Weighted_point_3                            Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Weighted_point_3>::value));

public:

  typedef Dimension_tag<3>  Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  typedef typename R_::Kernel_base::Weighted_point_3  Rep;
  typedef typename R_::Cartesian_const_iterator_3 Cartesian_const_iterator;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;

  typedef Point_3 Point;
  typedef FT Weight;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  Weighted_point_3() {}

  Weighted_point_3(const Origin& o)
    : Rep(typename R::Construct_weighted_point_3()(Return_base_tag(), o))
  {}

  Weighted_point_3(const Rep& p)
      : Rep(p) {}

  explicit
  Weighted_point_3(const Point_3& p)
    : Rep(typename R::Construct_weighted_point_3()(Return_base_tag(), p, 0))
  {}

  Weighted_point_3(const Point_3& p, const Weight& w)
    : Rep(typename R::Construct_weighted_point_3()(Return_base_tag(), p, w))
  {}

  Weighted_point_3(const FT& x, const FT& y, const FT& z)
    : Rep(typename R::Construct_weighted_point_3()(Return_base_tag(), x, y, z))
  {}

  typename cpp11::result_of<typename R::Construct_point_3(Weighted_point_3)>::type
  point() const
  {
    return typename R::Construct_point_3()(*this);
  }

  typename cpp11::result_of<typename R::Compute_weight_3(Weighted_point_3)>::type
  weight() const
  {
    return typename R::Compute_weight_3()(*this);
  }


  typename cpp11::result_of<typename R::Compute_x_3(Point_3)>::type
  x() const
  {
    return typename R::Compute_x_3()(point());
  }

  typename cpp11::result_of<typename R::Compute_y_3(Point_3)>::type
  y() const
  {
    return typename R::Compute_y_3()(point());
  }

  typename cpp11::result_of<typename R::Compute_z_3(Point_3)>::type
  z() const
  {
    return typename R::Compute_z_3()(point());
  }

  typename cpp11::result_of<typename R::Compute_hx_3(Point_3)>::type
  hx() const
  {
    return R().compute_hx_3_object()(point());
  }

  typename cpp11::result_of<typename R::Compute_hy_3(Point_3)>::type
  hy() const
  {
    return R().compute_hy_3_object()(point());
  }

  typename cpp11::result_of<typename R::Compute_hz_3(Point_3)>::type
  hz() const
  {
    return R().compute_hz_3_object()(point());
  }

  typename cpp11::result_of<typename R::Compute_hw_3(Point_3)>::type
  hw() const
  {
    return R().compute_hw_3_object()(point());
  }

  typename cpp11::result_of<typename R::Compute_x_3(Point_3)>::type
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

  typename cpp11::result_of<typename R::Compute_x_3(Point_3)>::type
  operator[](int i) const
  {
      return cartesian(i);
  }

  Cartesian_const_iterator cartesian_begin() const
  {
    return typename R::Construct_cartesian_const_iterator_3()(point());
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return typename R::Construct_cartesian_const_iterator_3()(point(),3);
  }

  int dimension() const
  {
      return 3;
  }

  Bbox_3 bbox() const
  {
    return R().construct_bbox_3_object()(point());
  }

  Weighted_point_3 transform(const Aff_transformation_3 &t) const
  {
    return Weighted_point_3(t.transform(point()),weight());
  }

};

template <class R>
inline
bool
operator==(const Origin& o, const Weighted_point_3<R>& p)
{ return p == o; }

template <class R>
inline
bool
operator!=(const Origin& o, const Weighted_point_3<R>& p)
{ return p != o; }

template <class R>
inline
bool
operator==(const Point_3<R>& bp, const Weighted_point_3<R>& p)
{ return bp == p.point(); }

template <class R>
inline
bool
operator!=(const Point_3<R>& bp, const Weighted_point_3<R>& p)
{ return bp != p.point(); }

template <class R>
inline
bool
operator==(const Weighted_point_3<R>& p, const Point_3<R>& bp)
{ return bp == p.point(); }

template <class R>
inline
bool
operator!=(const Weighted_point_3<R>& p, const Point_3<R>& bp)
{ return bp != p.point(); }

template <class R>
inline
bool
operator==(const Weighted_point_3<R>& p, const Weighted_point_3<R>& p2)
{ return p.point() == p2.point(); }

template <class R>
inline
bool
operator!=(const Weighted_point_3<R>& p, const Weighted_point_3<R>& p2)
{ return p.point() != p2.point(); }


template <class R>
inline
bool
operator<(const Weighted_point_3<R>& p, const Weighted_point_3<R>& q)
{ return p.point() < q.point(); }


template <class R >
std::ostream&
insert(std::ostream& os, const Weighted_point_3<R>& p,const Cartesian_tag&)
{
    switch(get_mode(os)) {
    case IO::ASCII :
        return os << p.point() << ' ' << p.weight();
    case IO::BINARY :
        write(os, p.x());
        write(os, p.y());
        write(os, p.z());
        write(os, p.weight());
        return os;
    default:
        return os << "Weighted_pointC3(" << p.x() << ", " << p.y()
                  << ", " << p.z() << ", " << p.weight() << ')';
    }
}

template <class R >
std::ostream&
insert(std::ostream& os, const Weighted_point_3<R>& p,const Homogeneous_tag&)
{
  switch(get_mode(os))
  {
    case IO::ASCII :
        return os << p.hx() << ' ' << p.hy() << ' ' << p.hz() << ' ' << p.hw() << ' ' << p.weight();
    case IO::BINARY :
        write(os, p.hx());
        write(os, p.hy());
        write(os, p.hz());
        write(os, p.hw());
        write(os, p.weight());
        return os;
    default:
        return os << "Weighted_pointH3(" << p.hx() << ", "
                                << p.hy() << ", "
                                << p.hz() << ", "
                                << p.hw() << ", "
                                << p.weight() << ')';
  }
}

template < class R >
std::ostream&
operator<<(std::ostream& os, const Weighted_point_3<R>& p)
{
  return insert(os, p, typename R::Kernel_tag() );
}


template <class R >
std::istream&
extract(std::istream& is, Weighted_point_3<R>& p, const Cartesian_tag&)
{
  typename R::FT x, y, z, weight;
    switch(get_mode(is)) {
    case IO::ASCII :
      is >> iformat(x) >> iformat(y) >> iformat(z) >> iformat(weight);
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        read(is, z);
        read(is, weight);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
      p = Weighted_point_3<R>(Point_3<R>(x, y, z), weight);
    return is;
}


template <class R >
std::istream&
extract(std::istream& is, Weighted_point_3<R>& p, const Homogeneous_tag&)
{
  typename R::RT hx, hy, hz, hw;
  typename R::FT weight;
  switch(get_mode(is))
  {
    case IO::ASCII :
      is >> hx >> hy >> hz >> hw >> weight;
        break;
    case IO::BINARY :
        read(is, hx);
        read(is, hy);
        read(is, hz);
        read(is, hw);
        read(is, weight);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  if (is)
    p = Weighted_point_3<R>(Point_3<R>(hx, hy, hz, hw), weight);
  return is;
}

template < class R >
std::istream&
operator>>(std::istream& is, Weighted_point_3<R>& p)
{
  return extract(is, p, typename R::Kernel_tag() );
}

} //namespace CGAL

#endif // CGAL_WEIGHTED_POINT_3_H
