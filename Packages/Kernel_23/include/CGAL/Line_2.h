// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_LINE_2_H
#define CGAL_LINE_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Line_2 : public R_::Kernel_base::Line_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Segment_2             Segment_2;
  typedef typename R_::Ray_2                 Ray_2;
  typedef typename R_::Vector_2              Vector_2;
  typedef typename R_::Direction_2           Direction_2;
  typedef typename R_::Kernel_base::Line_2  RLine_2;
public:
  typedef  R_   R;

  Line_2() {}

  Line_2(const Point_2 &p, const Point_2 &q)
    : RLine_2(p,q) {}

  Line_2(const RT &a, const RT &b, const RT &c)
    : RLine_2(a,b,c) {}

  Line_2(const RLine_2& l)  // conversion impl -> interface class
    : RLine_2(l) {}

  Line_2(const Segment_2& s)
    : RLine_2(s) {}

  Line_2(const Ray_2& r)
    : RLine_2(r) {}

  Line_2(const Point_2 &p, const Direction_2 &d)
    : RLine_2(p,d) {}

  Line_2(const Point_2 &p, const Vector_2 &v)
    : RLine_2(p,v) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_LINE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Line_2<R> &l)
{
  typedef typename  R::Kernel_base::Line_2  RLine_2;
  return os << static_cast<const RLine_2&>(l);
}
#endif // CGAL_NO_OSTREAM_INSERT_LINE_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINE_2
template < class R >
std::istream &
operator>>(std::istream &is, Line_2<R> &p)
{
  typedef typename  R::Kernel_base::Line_2  RLine_2;
  return is >> static_cast<RLine_2&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINE_2

CGAL_END_NAMESPACE

#endif  // CGAL_LINE_2_H
