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

#ifndef CGAL_CARTESIAN_SEGMENT_2_H
#define CGAL_CARTESIAN_SEGMENT_2_H

#include <CGAL/Twotuple.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>

CGAL_BEGIN_NAMESPACE
template < class R_ >
class SegmentC2
{
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Twotuple<Point_2>                        Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                     R;

  SegmentC2() 
  {}

  SegmentC2(const Point_2 &sp, const Point_2 &ep)
    : base(sp, ep) 
  {}

  const Point_2 &   
  source() const
  {
      return get(base).e0;
  }
  
  const Point_2 &   
  target() const
  {
      return get(base).e1;
  }

  Segment_2 transform(const Aff_transformation_2 &t) const
  {
    return Segment_2(t.transform(source()), t.transform(target()));
  }

};

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENTC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const SegmentC2<R> &s)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << s.source() << ' ' << s.target();
    case IO::BINARY :
        return os << s.source() << s.target();
    default:
        return os << "SegmentC2(" << s.source() <<  ", " << s.target() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENTC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENTC2
template < class R >
std::istream &
operator>>(std::istream &is, SegmentC2<R> &s)
{
    typename R::Point_2 p, q;

    is >> p >> q;

    if (is)
	s = SegmentC2<R>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENTC2






CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SEGMENT_2_H
