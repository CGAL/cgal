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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_PREDICATES_ON_POINTS_3_H
#define CGAL_CARTESIAN_PREDICATES_ON_POINTS_3_H

#include <CGAL/predicates/kernel_ftC3.h>

namespace CGAL {

template < class K >
inline
bool
equal_xy(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().equal_xy_3_object()(p, q);
}

template < class K >
inline
bool
equal_xyz(const PointC3<K> &p, const PointC3<K> &q)
{
  return p.x() == q.x() && p.y() == q.y() && p.z() == q.z();
}

template < class K >
inline
Comparison_result
compare_xy(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().compare_xy_3_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_lexicographically_xy(const PointC3<K> &p, const PointC3<K> &q)
{
  return K().compare_xy_3_object()(p, q);
}

template < class K >
inline
bool
lexicographically_xy_smaller_or_equal(const PointC3<K> &p, 
				      const PointC3<K> &q)
{ 
  return compare_lexicographically_xy(p, q) != LARGER;
}

template < class K >
inline
bool
lexicographically_xy_smaller(const PointC3<K> &p, 
			     const PointC3<K> &q)
{ 
  return K().less_xy_3_object()(p, q);
}

template < class K >
inline
bool
strict_dominance(const PointC3<K> &p,
		 const PointC3<K> &q)
{
  return strict_dominanceC3(p.x(), p.y(), p.z(),
			    q.x(), q.y(), q.z());
}

template < class K >
inline
bool
dominance(const PointC3<K> &p,
	  const PointC3<K> &q)
{
  return dominanceC3(p.x(), p.y(), p.z(),
		     q.x(), q.y(), q.z());
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_PREDICATES_ON_POINTS_3_H
