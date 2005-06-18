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

#ifndef CGAL_CARTESIAN_RAY_2_H
#define CGAL_CARTESIAN_RAY_2_H

#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class RayC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Twotuple<Point_2>                        Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  RayC2() 
  {}

  RayC2(const Point_2 &sp, const Point_2 &secondp)
    : base(sp, secondp) 
  {}


  const Point_2&
  source() const
  {
    return get(base).e0;
  }

  const Point_2 &
  second_point() const
  {
    return get(base).e1;
  }
  
  Ray_2 
  transform(const Aff_transformation_2 &t) const
  {
    return Ray_2(t.transform(source()), t.transform(second_point()));
  }
};


#ifndef CGAL_NO_OSTREAM_INSERT_RAYC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const RayC2<R> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r.source() << ' ' << r.second_point();
    case IO::BINARY :
        return os << r.source() << r.second_point();
    default:
        return os << "RayC2(" << r.source() <<  ", " << r.second_point() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_RAYC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAYC2
template < class R >
std::istream &
operator>>(std::istream &is, RayC2<R> &r)
{
    typename R::Point_2 p, q;

    is >> p >> q;

    if (is)
	r = RayC2<R>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAYC2


CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_RAY_2_H






