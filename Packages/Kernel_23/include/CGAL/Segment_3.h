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
// Author(s)     : Andreas Fabri, Stefan Schirra

#ifndef CGAL_SEGMENT_3_H
#define CGAL_SEGMENT_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Segment_3 : public R_::Kernel_base::Segment_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Kernel_base::Segment_3  RSegment_3;
public:
  typedef          R_                       R;

  Segment_3() {}

  Segment_3(const Point_3& sp, const Point_3& ep)
    : RSegment_3(sp,ep) {}

  Segment_3(const RSegment_3& s)
      : RSegment_3(s) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENT_3
template < class R>
std::ostream&
operator<<(std::ostream& os, const Segment_3<R>& s)
{
  typedef typename  R::Kernel_base::Segment_3  RSegment_3;
  return os << static_cast<const RSegment_3&>(s);
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENT_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENT_3
template < class R>
std::istream&
operator>>(std::istream& is, Segment_3<R>& s)
{
  typedef typename  R::Kernel_base::Segment_3  RSegment_3;
  return is >> static_cast<RSegment_3&>(s);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENT_3

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_3_H
