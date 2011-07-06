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
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_POINT_2_LINE_2_INTERSECTION_H
#define CGAL_POINT_2_LINE_2_INTERSECTION_H

#include <CGAL/Line_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Object.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {

namespace internal {

template <class K>
inline bool
do_intersect(const typename K::Point_2 &pt, 
	     const typename K::Line_2 &line,
	     const K&)
{
    return line.has_on(pt);
}

template <class K>
inline bool
do_intersect(const typename K::Line_2 &line,
	     const typename K::Point_2 &pt, 
	     const K&)
{
    return line.has_on(pt);
}

template <class K>
typename CGAL::Intersection_traits_2
<K, typename K::Point_2, typename K::Line_2>::result_type
intersection(const typename K::Point_2 &pt, 
	     const typename K::Line_2 &line,
	     const K& k)
{
  typedef typename CGAL::Intersection_traits_2
    <K, typename K::Point_2, typename K::Line_2>::result_type result_type;

    if (do_intersect(pt,line, k)) {
        return result_type(pt);
    }
    return result_type();
}

template <class K>
typename CGAL::Intersection_traits_2
<K, typename K::Point_2, typename K::Line_2>::result_type
intersection(const typename K::Line_2 &line,
	     const typename K::Point_2 &pt, 
	     const K& k)
{
  return internal::intersection(pt, line, k);
}

} // namespace internal

} //namespace CGAL

#endif
