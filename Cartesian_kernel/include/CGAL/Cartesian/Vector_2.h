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
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_VECTOR_2_H
#define CGAL_CARTESIAN_VECTOR_2_H

#include <CGAL/Origin.h>
#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class VectorC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Twotuple<FT>	                           Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

  static FT one;

public:
  typedef R_                                     R;

  VectorC2() {}
 
  VectorC2(const FT &x, const FT &y)
    : base(x, y) {}

  VectorC2(const FT &hx, const FT &hy, const FT &hw)
  {
    if (hw != FT(1))
      base = Rep(hx/hw, hy/hw);
    else
      base = Rep(hx, hy);
  }

  const FT & x() const
  {
      return get(base).e0;
  }
  const FT & y() const
  {
      return get(base).e1;
  }

  const FT & hx() const
  {
      return x();
  }
  const FT & hy() const
  {
      return y();
  }
  const FT& hw() const
  {
    return one;
  }
  


  Vector_2 transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }
};

template <class R >
typename R::FT VectorC2<R>::one = 1;

template < class R >
CGAL_KERNEL_INLINE
bool 
operator==(const VectorC2<R> &v, const VectorC2<R> &w)
{
  return w.x() == v.x() && w.y() == v.y();
}

template < class R >
inline
bool 
operator!=(const VectorC2<R> &v, const VectorC2<R> &w)
{
  return !(v == w);
}

template < class R >
inline
bool
operator==(const VectorC2<R> &v, const Null_vector &)
{
  return CGAL_NTS is_zero(v.x()) && CGAL_NTS is_zero(v.y());
}

template < class R >
inline
bool
operator==(const Null_vector &n, const VectorC2<R> &v)
{
  return v == n;
}

template < class R >
inline
bool
operator!=(const VectorC2<R> &v, const Null_vector &n)
{
  return !(v == n);
}

template < class R >
inline
bool
operator!=(const Null_vector &n, const VectorC2<R> &v)
{
  return !(v == n);
}



#ifndef CGAL_NO_OSTREAM_INSERT_VECTORC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const VectorC2<R> &v)
{
    switch(os.iword(IO::mode)) {
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
#endif // CGAL_NO_OSTREAM_INSERT_VECTORC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTORC2
template < class R >
std::istream &
operator>>(std::istream &is, VectorC2<R> &p)
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
	p = VectorC2<R>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTORC2

CGAL_END_NAMESPACE
#endif // CGAL_CARTESIAN_VECTOR_2_H
