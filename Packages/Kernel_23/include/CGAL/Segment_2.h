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

#ifndef CGAL_SEGMENT_2_H
#define CGAL_SEGMENT_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Segment_2 : public R_::Kernel_base::Segment_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Kernel_base::Segment_2  RSegment_2;
public:
  typedef  R_                               R;

  Segment_2() {}

  Segment_2(const Point_2 &sp, const Point_2 &ep)
    :  RSegment_2(sp,ep) {}

  // conversion from implementation class object to interface class object
  Segment_2(const RSegment_2& s)
    : RSegment_2(s) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_SEGMENT_2
template < class R>
std::ostream &
operator<<(std::ostream &os, const Segment_2<R> &s)
{
  typedef typename  R::Kernel_base::Segment_2  RSegment_2;
  return os << static_cast<const RSegment_2&>(s);
}
#endif // CGAL_NO_OSTREAM_INSERT_SEGMENT_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_SEGMENT_2
template < class R>
std::istream &
operator>>(std::istream &is, Segment_2<R> &s)
{
  typedef typename  R::Kernel_base::Segment_2  RSegment_2;
  return is >> static_cast<RSegment_2&>(s);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SEGMENT_2

CGAL_END_NAMESPACE

#endif //  CGAL_SEGMENT_2_H
