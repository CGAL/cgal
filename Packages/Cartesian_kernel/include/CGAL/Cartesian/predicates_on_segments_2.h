// Copyright (c) 2000  Utrecht University (The Netherlands),
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
// Author(s)     : Susan Hert

#ifndef CGAL_CARTESIAN_PREDICATES_ON_SEGMENTS_2_H
#define CGAL_CARTESIAN_PREDICATES_ON_SEGMENTS_2_H

#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/Cartesian/Segment_2.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
Comparison_result
compare_slopes(const SegmentC2<K> &s1, const SegmentC2<K> &s2)
{
  return K().compare_slope_2_object()(s1, s2);
}

template < class K >
inline
Comparison_result
compare_y_at_x(const PointC2<K> &p, const SegmentC2<K> &s)
{
  return K().compare_y_at_x_2_object()(p, s);
}

template < class K >
inline
Comparison_result
compare_y_at_x(const PointC2<K> &p,
               const SegmentC2<K> &s1,
               const SegmentC2<K> &s2)
{
  return K().compare_y_at_x_2_object()(p, s1, s2);
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_SEGMENTS_2_H
