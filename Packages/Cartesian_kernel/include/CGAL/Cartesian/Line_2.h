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

#ifndef CGAL_CARTESIAN_LINE_2_H
#define CGAL_CARTESIAN_LINE_2_H

#include <CGAL/Threetuple.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class LineC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Line_2               Line_2;


  typedef Threetuple<FT>	                   Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  LineC2() {}

  LineC2(const FT &a, const FT &b, const FT &c)
    : base(a, b, c) {}
  
  bool            operator==(const LineC2 &l) const;
  bool            operator!=(const LineC2 &l) const;

  const FT & a() const
  {
      return get(base).e0;
  }
  const FT & b() const
  {
      return get(base).e1;
  }
  const FT & c() const
  {
      return get(base).e2;
  }



};

template < class R >
CGAL_KERNEL_INLINE
bool
LineC2<R>::operator==(const LineC2<R> &l) const
{
  if (CGAL::identical(base, l.base))
      return true;
  return equal_lineC2(a(), b(), c(), l.a(), l.b(), l.c());
}

template < class R >
inline
bool
LineC2<R>::operator!=(const LineC2<R> &l) const
{
  return ! (*this == l);
}






#ifndef CGAL_NO_OSTREAM_INSERT_LINEC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const LineC2<R> &l)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << l.a() << ' ' << l.b() << ' ' << l.c();
    case IO::BINARY :
        write(os, l.a());
        write(os, l.b());
        write(os, l.c());
        return os;
    default:
        return os << "LineC2(" << l.a() 
		  << ", " << l.b() << ", " << l.c() <<')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_LINEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINEC2
template < class R >
std::istream &
operator>>(std::istream &is, LineC2<R> &l)
{
    typename R::FT a, b, c;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> a >> b >> c;
        break;
    case IO::BINARY :
        read(is, a);
        read(is, b);
        read(is, c);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	l = LineC2<R>(a, b, c);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINEC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_LINE_2_H
