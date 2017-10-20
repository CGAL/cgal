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

#ifndef CGAL_CARTESIAN_PREDICATES_ON_DIRECTIONS_2_H
#define CGAL_CARTESIAN_PREDICATES_ON_DIRECTIONS_2_H

namespace CGAL {

template < class K >
inline
typename K::Boolean
equal_direction(const DirectionC2<K> &d1,
                const DirectionC2<K> &d2)
{
  return equal_directionC2(d1.dx(), d1.dy(), d2.dx(), d2.dy());
}

template < class K >
inline
typename K::Comparison_result
compare_angle_with_x_axis(const DirectionC2<K> &d1,
                          const DirectionC2<K> &d2)
{
  return K().compare_angle_with_x_axis_2_object()(d1, d2);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_PREDICATES_ON_DIRECTIONS_2_H
