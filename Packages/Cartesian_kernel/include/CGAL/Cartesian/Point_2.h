// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
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
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_POINT_2_H
#define CGAL_CARTESIAN_POINT_2_H

#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Twotuple.h>
#include <CGAL/Kernel/Cartesian_coordinate_iterator_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointC2
  : public R_::template Handle<Twotuple<typename R_::FT> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Twotuple<FT>	                           rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef Cartesian_coordinate_iterator_2<R_> Cartesian_const_iterator;

  typedef R_                                     R;

  PointC2() {}

  PointC2(const Origin &)
    : base(rep(FT(0), FT(0))) {}

  PointC2(const FT &x, const FT &y)
    : base(rep(x, y)) {}

  PointC2(const FT &hx, const FT &hy, const FT &hw)
  {
    if (hw != FT(1))
      initialize_with( rep(hx/hw, hy/hw) );
    else
      initialize_with( rep(hx, hy) );
  }

  PointC2(const Vector_2 &v)
    : base(v) {}

  const FT& x() const
  {
      return Ptr()->e0;
  }
  const FT& y() const
  {
      return Ptr()->e1;
  }

  const FT& hx() const
  {
      return x();
  }
  const FT& hy() const
  {
      return y();
  }
  FT hw() const
  {
      return FT(1);
  }

  const FT& cartesian(int i) const;
  FT homogeneous(int i) const;
  const FT& operator[](int i) const
  {
      return cartesian(i);
  }


  Cartesian_const_iterator cartesian_begin() const 
  {
    return Cartesian_const_iterator(static_cast<const Point_2* >(this),0);
  }

  Cartesian_const_iterator cartesian_end() const 
  {
    return Cartesian_const_iterator(static_cast<const Point_2* >(this), 2);
  }

  int dimension() const
  {
      return 2;
  }

  bool operator==(const PointC2 &p) const
  {
      if (identical(p))
	  return true;
      return equal_xy(*this, p);
  }
  bool operator!=(const PointC2 &p) const
  {
      return !(*this == p);
  }

  Bbox_2 bbox() const;

  Point_2 transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }
};

template < class R >
CGAL_KERNEL_INLINE
const typename PointC2<R>::FT &
PointC2<R>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i == 0) || (i == 1) );
  return (i == 0) ? x() : y();
}

template < class R >
CGAL_KERNEL_INLINE
typename PointC2<R>::FT
PointC2<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  if (i<2)
    return cartesian(i);
  return FT(1);
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_2
PointC2<R>::bbox() const
{
  std::pair<double,double> xp = CGAL::to_interval(x());
  std::pair<double,double> yp = CGAL::to_interval(y());
  return Bbox_2(xp.first, yp.first,  xp.second, yp.second);
}

#ifndef CGAL_NO_OSTREAM_INSERT_POINTC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const PointC2<R> &p)
{
    switch(os.iword(IO::mode)) {
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
#endif // CGAL_NO_OSTREAM_INSERT_POINTC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINTC2
template < class R >
std::istream &
operator>>(std::istream &is, PointC2<R> &p)
{
    typename R::FT x, y;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> x >> y;
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	p = PointC2<R>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINTC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_2_H
