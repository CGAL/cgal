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

#ifndef CGAL_CARTESIAN_TRIANGLE_2_H
#define CGAL_CARTESIAN_TRIANGLE_2_H

#include <CGAL/Threetuple.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class TriangleC2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Triangle_2           Triangle_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Threetuple<Point_2>	                   Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                    R;

  TriangleC2() {}

  TriangleC2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
    : base(p, q, r) {}


  const Point_2 &
  vertex(int i) const
  {
    if (i>2) i = i%3;
    else if (i<0) i = (i%3) + 3;
    return (i==0) ? get(base).e0 :
      (i==1) ? get(base).e1 :
      get(base).e2;
  }
  
  const Point_2 &
  operator[](int i) const
  {
    return vertex(i);
  }

  Triangle_2           
  transform(const Aff_transformation_2 &t) const
  {
    return Triangle_2(t.transform(vertex(0)),
		      t.transform(vertex(1)),
		      t.transform(vertex(2)));
  }
  
};

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLEC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const TriangleC2<R> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0] << t[1]  << t[2];
    default:
        return os<< "TriangleC2(" << t[0] << ", " 
		 << t[1] << ", " << t[2] <<")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLEC2
template < class R >
std::istream &
operator>>(std::istream &is, TriangleC2<R> &t)
{
    typename R::Point_2 p, q, r;

    is >> p >> q >> r;

    if (is)
	t = TriangleC2<R>(p, q, r);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLEC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TRIANGLE_2_H
