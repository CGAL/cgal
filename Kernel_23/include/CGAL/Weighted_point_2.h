// Copyright (c) 2016
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
// Author(s)     : Mariette Yvinec, Sylvain Pion

#ifndef CGAL_WEIGHTED_POINT_2_H
#define CGAL_WEIGHTED_POINT_2_H

#include <CGAL/Origin.h>
#include <CGAL/representation_tags.h>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R_>
class Weighted_point_2 : public R_::Kernel_base::Weighted_point_2
{
  typedef typename R_::FT                             FT;
  typedef typename R_::FT                             RT;

  typedef Weighted_point_2<R_>                        Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Weighted_point_2>::value));

public:
  typedef Dimension_tag<2>                            Ambient_dimension;
  typedef Dimension_tag<0>                            Feature_dimension;

  typedef typename R_::Kernel_base::Weighted_point_2  Rep;
  typedef typename R_::Cartesian_const_iterator_2     Cartesian_const_iterator;
  typedef typename R_::Point_2                        Point_2;
  typedef typename R_::Aff_transformation_2           Aff_transformation_2;

  typedef Point_2                                     Point;
  typedef FT                                          Weight;
  typedef          R_                                 R;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  Weighted_point_2() {}

  Weighted_point_2(const Origin& o)
    : Rep(typename R::Construct_weighted_point_2()(Return_base_tag(), o))
  {}

  Weighted_point_2(const Rep& p)
      : Rep(p) {}

  explicit
  Weighted_point_2(const Point_2& p)
    : Rep(typename R::Construct_weighted_point_2()(Return_base_tag(), p, 0))
  {}

  Weighted_point_2(const Point_2& p, const Weight& w)
    : Rep(typename R::Construct_weighted_point_2()(Return_base_tag(), p, w))
  {}

  Weighted_point_2(const FT& x, const FT& y)
    : Rep(typename R::Construct_weighted_point_2()(Return_base_tag(), x, y))
  {}

  typename cpp11::result_of<typename R::Construct_point_2( Weighted_point_2)>::type
  point() const
  {
    return typename R::Construct_point_2()(*this);
  }

  typename cpp11::result_of<typename R::Compute_weight_2( Weighted_point_2)>::type
  weight() const
  {
    return typename R::Compute_weight_2()(*this);
  }


  typename cpp11::result_of<typename R::Compute_x_2( Point_2)>::type
  x() const
  {
    return typename R::Compute_x_2()(point());
  }

  typename cpp11::result_of<typename R::Compute_y_2( Point_2)>::type
  y() const
  {
    return typename R::Compute_y_2()(point());
  }

  typename cpp11::result_of<typename R::Compute_hx_2( Point_2)>::type
  hx() const
  {
    return R().compute_hx_2_object()(point());
  }

  typename cpp11::result_of<typename R::Compute_hy_2( Point_2)>::type
  hy() const
  {
    return R().compute_hy_2_object()(point());
  }

  typename cpp11::result_of<typename R::Compute_hw_2( Point_2)>::type
  hw() const
  {
    return R().compute_hw_2_object()(point());
  }

  typename cpp11::result_of<typename R::Compute_x_2( Point_2)>::type
  cartesian(int i) const
  {
    CGAL_kernel_precondition( (i == 0) || (i == 1) );
    if (i==0) return x();
    return y();
  }

  RT
  homogeneous(int i) const
  {
    CGAL_kernel_precondition( (i >= 0) || (i <= 2) );
    if (i==0) return hx();
    if (i==1) return hy();
    return hw();
  }

  typename cpp11::result_of<typename R::Compute_x_2(Point_2)>::type
  operator[](int i) const
  {
      return cartesian(i);
  }

  Cartesian_const_iterator cartesian_begin() const
  {
    return typename R::Construct_cartesian_const_iterator_2()(point());
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return typename R::Construct_cartesian_const_iterator_2()(point(),3);
  }

  int dimension() const
  {
      return 2;
  }

  Bbox_2 bbox() const
  {
    return R().construct_bbox_2_object()(point());
  }

  Weighted_point_2 transform(const Aff_transformation_2 &t) const
  {
    return Weighted_point_2(t.transform(point(),weight()));
  }

};

template <class R>
inline
bool
operator==(const Origin& o, const Weighted_point_2<R>& p)
{ return p == o; }

template <class R>
inline
bool
operator!=(const Origin& o, const Weighted_point_2<R>& p)
{ return p != o; }


template <class R>
inline
bool
operator==(const Point_2<R>& bp, const Weighted_point_2<R>& p)
{ return bp == p.point(); }

template <class R>
inline
bool
operator!=(const Point_2<R>& bp, const Weighted_point_2<R>& p)
{ return bp != p.point(); }

template <class R>
inline
bool
operator==(const Weighted_point_2<R>& p, const Point_2<R>& bp)
{ return bp == p.point(); }

template <class R>
inline
bool
operator!=(const Weighted_point_2<R>& p, const Point_2<R>& bp)
{ return bp != p.point(); }

template <class R>
inline
bool
operator==(const Weighted_point_2<R>& p, const Weighted_point_2<R>& p2)
{ return p.point() == p2.point(); }

template <class R>
inline
bool
operator!=(const Weighted_point_2<R>& p, const Weighted_point_2<R>& p2)
{ return p.point() != p2.point(); }


template <class R>
inline
bool
operator<(const Weighted_point_2<R>& p, const Weighted_point_2<R>& q)
{ return p.point() < q.point(); }


template <class R >
std::ostream&
insert(std::ostream& os, const Weighted_point_2<R>& p,const Cartesian_tag&)
{
    switch(get_mode(os)) {
    case IO::ASCII :
        return os << p.point() << ' ' << p.weight();
    case IO::BINARY :
        write(os, p.x());
        write(os, p.y());
        write(os, p.weight());
        return os;
    default:
      return os << "Weighted_pointC2(" << p.x() << ", " << p.y() << ", " << p.weight()<< ')';
    }
}

template <class R >
std::ostream&
insert(std::ostream& os, const Weighted_point_2<R>& p,const Homogeneous_tag&)
{
  switch(get_mode(os))
  {
    case IO::ASCII :
      return os << p.point() << ' ' << p.weight();
    case IO::BINARY :
      os << p.point();
      write(os, p.weight());
      return os;
    default:
      return os << "Weighted_pointH2("
                << p.hx() << ", "
                << p.hy() << ", "
                << p.hw() << ", "
                << p.weight() << ')';
  }
}

template < class R >
std::ostream&
operator<<(std::ostream& os, const Weighted_point_2<R>& p)
{
  return insert(os, p, typename R::Kernel_tag() );
}


template <class R >
std::istream&
extract(std::istream& is, Weighted_point_2<R>& p, const Cartesian_tag&)
{
  typename R::FT x, y, weight;
    switch(get_mode(is)) {
    case IO::ASCII :
      is >> iformat(x) >> iformat(y) >> iformat(weight);
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        read(is, weight);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
      p = Weighted_point_2<R>(typename R::Point_2(x, y),weight);
    return is;
}


template <class R >
std::istream&
extract(std::istream& is, Weighted_point_2<R>& p, const Homogeneous_tag&)
{
  typename R::RT hx, hy, hw;
  typename R::FT weight;
  switch(get_mode(is))
  {
    case IO::ASCII :
      is >> hx >> hy >> hw >> weight;
        break;
    case IO::BINARY :
        read(is, hx);
        read(is, hy);
        read(is, hw);
        read(is, weight);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  if (is)
    p = Weighted_point_2<R>(typename R::Point_2(hx, hy, hw),weight);
  return is;
}

template < class R >
std::istream&
operator>>(std::istream& is, Weighted_point_2<R>& p)
{
  return extract(is, p, typename R::Kernel_tag() );
}

} //namespace CGAL

#endif // CGAL_WEIGHTED_POINT_2_H
