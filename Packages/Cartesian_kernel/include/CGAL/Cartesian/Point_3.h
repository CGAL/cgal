// Copyright (c) 2000  Utrecht University (The Netherlands),
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
// Author(s)     : Andreas Fabri and Herve Bronnimann

#ifndef CGAL_CARTESIAN_POINT_3_H
#define CGAL_CARTESIAN_POINT_3_H

#include <CGAL/Threetuple.h>
#include <CGAL/Origin.h>
#include <CGAL/Bbox_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointC3
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef Threetuple<FT>                           Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef const FT* Cartesian_const_iterator;
  typedef R_                                R;

  PointC3() {}

  PointC3(const Origin &)
    : base(FT(0), FT(0), FT(0)) {}

  PointC3(const FT &x, const FT &y, const FT &z)
    : base(x, y, z) {}

  PointC3(const FT &x, const FT &y, const FT &z, const FT &w)
  {
    if (w != FT(1))
      base = Rep(x/w, y/w, z/w);
    else
      base = Rep(x, y, z);
  }

/*
  bool operator==(const PointC3 &p) const
  {
      if (CGAL::identical(base, p.base))
	  return true;
      return x_equal(*this, p) && y_equal(*this, p) && z_equal(*this, p);
  }
  bool operator!=(const PointC3 &p) const
  {
      return !(*this == p);
  }
*/

  const FT & x() const
  {
      return get(base).e0;
  }
  const FT & y() const
  {
      return get(base).e1;
  }
  const FT & z() const
  {
      return get(base).e2;
  }

  const FT & hx() const
  {
      return x();
  }
  const FT & hy() const
  {
      return y();
  }
  const FT & hz() const
  {
      return z();
  }
  FT hw() const
  {
      return FT(1);
  }

  const FT & cartesian(int i) const;
  const FT & operator[](int i) const;
  FT homogeneous(int i) const;

  Cartesian_const_iterator cartesian_begin() const 
  {
    return & get(base).e0; 
  }

  Cartesian_const_iterator cartesian_end() const 
  {
    const FT* ptr = & get(base).e2;
    ptr++;
    return ptr;
  }

  int dimension() const
  {
      return 3;
  }
  Bbox_3 bbox() const;

  Point_3 transform(const Aff_transformation_3 &t) const
  {
    return t.transform(*this);
  }
};

template < class R >
inline
const typename PointC3<R>::FT &
PointC3<R>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  return *(&(get(base).e0)+i);
}

template < class R >
inline
const typename PointC3<R>::FT &
PointC3<R>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
inline
typename PointC3<R>::FT
PointC3<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition(i>=0 && i<=3);
  if (i<3) return cartesian(i);
  return FT(1);
}

template < class R >
Bbox_3
PointC3<R>::bbox() const
{
  std::pair<double,double> xp = CGAL_NTS to_interval(x());
  std::pair<double,double> yp = CGAL_NTS to_interval(y());
  std::pair<double,double> zp = CGAL_NTS to_interval(z());
  return Bbox_3(xp.first, yp.first, zp.first, xp.second, yp.second, zp.second);
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_POINTC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const PointC3<R> &p)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << p.x() << ' ' << p.y()  << ' ' << p.z();
    case IO::BINARY :
        write(os, p.x());
        write(os, p.y());
        write(os, p.z());
        return os;
    default:
        os << "PointC3(" << p.x() << ", " << p.y()  << ", " << p.z() << ")";
        return os;
    }
}
#endif // CGAL_CARTESIAN_NO_OSTREAM_INSERT_POINTC3

#ifndef CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_POINTC3
template < class R >
std::istream &
operator>>(std::istream &is, PointC3<R> &p)
{
    typename R::FT x, y, z;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> x >> y >> z;
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
	p = PointC3<R>(x, y, z);
    return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_POINTC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_3_H
