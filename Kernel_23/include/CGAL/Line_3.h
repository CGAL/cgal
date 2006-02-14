// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// Author(s)     : Andreas Fabri
//                 Stefan Schirra

#ifndef CGAL_LINE_3_H
#define CGAL_LINE_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Line_3 : public R_::Kernel_base::Line_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Kernel_base::Line_3   RLine_3;
public:
  typedef          R_                       R;

  Line_3() {}

  Line_3(const Point_3 & p, const Point_3 & q)
      : RLine_3(p,q) {}

  // conversion impl -> interface class
  Line_3(const RLine_3& l)
      : RLine_3(l) {}

  Line_3(const Segment_3 & s)
      : RLine_3( s ) {}

  Line_3(const Ray_3 & r)
      : RLine_3( r ) {}

  Line_3(const Point_3 & p, const Direction_3 & d)
      : RLine_3( p, d ) {}

  Line_3(const Point_3 & p, const Vector_3 & v)
      : RLine_3( p, v ) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_LINE_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Line_3<R>& l)
{
  typedef typename  R::Kernel_base::Line_3  RLine_3;
  return os << static_cast<const RLine_3&>(l);
}
#endif // CGAL_NO_OSTREAM_INSERT_LINE_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINE_3
template < class R >
std::istream&
operator>>(std::istream & is, Line_3<R> & p)
{
  typedef typename  R::Kernel_base::Line_3  RLine_3;
  is >> static_cast<RLine_3&>(p);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINE_3

CGAL_END_NAMESPACE

#endif // CGAL_LINE_3_H
