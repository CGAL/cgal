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
// 
//
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_LINE_CONSTRUCTIONS_2_H
#define CGAL_CARTESIAN_LINE_CONSTRUCTIONS_2_H

#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/Cartesian/Line_2.h>

#include <CGAL/Cartesian/Direction_2.h>

namespace CGAL {

template < class K >
inline
LineC2<K>
line_from_point_direction(const PointC2<K> &p,
                          const DirectionC2<K> &d)
{
  return K().construct_line_2_object()(p, d);
}

template < class K >
inline
LineC2<K>
perpendicular_through_point(const LineC2<K> &l,
                            const PointC2<K> &p)
{
  typename K::FT a, b, c;
  perpendicular_through_pointC2(l.a(), l.b(), p.x(), p.y(), a, b, c);
  return LineC2<K>(a, b, c);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_LINE_CONSTRUCTIONS_2_H
