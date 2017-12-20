// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_PREDICATES_ON_POINTS_2_H
#define CGAL_CARTESIAN_PREDICATES_ON_POINTS_2_H

#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/predicates/kernel_ftC2.h>

namespace CGAL {

template < class K >
inline
bool
equal_xy(const PointC2<K> &p, const PointC2<K> &q)
{
  return CGAL_AND( p.x() == q.x() , p.y() == q.y() );
}

#if 0
// Unused, undocumented, un-functorized.
template < class K >
inline
Comparison_result
compare_deltax_deltay(const PointC2<K>& p,
                      const PointC2<K>& q,
                      const PointC2<K>& r,
                      const PointC2<K>& s)
{
  return compare_deltax_deltayC2(p.x(), q.x(), r.y(), s.y());
}
#endif

template < class K >
inline
Comparison_result
compare_lexicographically_yx(const PointC2<K> &p,
                             const PointC2<K> &q)
{
  return compare_lexicographically_xyC2(p.y(), p.x(), q.y(), q.x());
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_PREDICATES_ON_POINTS_2_H
