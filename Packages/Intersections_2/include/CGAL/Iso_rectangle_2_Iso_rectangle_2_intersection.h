
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
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_ISO_RECTANGLE_2_ISO_RECTANGLE_2_INTERSECTION_H
#define CGAL_ISO_RECTANGLE_2_ISO_RECTANGLE_2_INTERSECTION_H

#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
Object
intersection(
    const typename CGAL_WRAP(K)::Iso_rectangle_2 &irect1,
    const typename CGAL_WRAP(K)::Iso_rectangle_2 &irect2,
    const K&)
{
    typename K::Construct_point_2 construct_point_2;
    typename K::Construct_object_2 construct_object;
    typename K::Construct_iso_rectangle_2 construct_iso_rectangle_2;
    const typename K::Point_2 &min1 = irect1.min();
    const typename K::Point_2 &min2 = irect2.min();
    const typename K::Point_2 &max1 = irect1.max();
    const typename K::Point_2 &max2 = irect2.max();
    typename K::FT minx, miny, maxx, maxy;
    typename K::Point_2 newmin;
    typename K::Point_2 newmax;
    minx = (min1.x() >= min2.x()) ? min1.x() : min2.x();
    maxx = (max1.x() <= max2.x()) ? max1.x() : max2.x();
    if (maxx < minx)
        return Object();
    miny = (min1.y() >= min2.y()) ? min1.y() : min2.y();
    maxy = (max1.y() <= max2.y()) ? max1.y() : max2.y();
    if (maxy < miny)
        return Object(); 
    if (K::FT_denominator(minx) == K::FT_denominator(miny)) {
        newmin = construct_point_2(K::FT_numerator(minx), K::FT_numerator(miny),
				   K::FT_denominator(minx));
    } else {
        newmin = construct_point_2(K::FT_numerator(minx)*K::FT_denominator(miny),
				   K::FT_numerator(miny)*K::FT_denominator(minx),
				   K::FT_denominator(minx) * K::FT_denominator(miny));
    }
    if (K::FT_denominator(maxx) == K::FT_denominator(maxy)) {
        newmax = construct_point_2(K::FT_numerator(maxx), K::FT_numerator(maxy),
				   K::FT_denominator(maxx));
    } else {
        newmax = construct_point_2(K::FT_numerator(maxx)*K::FT_denominator(maxy),
				   K::FT_numerator(maxy)*K::FT_denominator(maxx),
				   K::FT_denominator(maxx) * K::FT_denominator(maxy));
    }
    return construct_object(construct_iso_rectangle_2(newmin, newmax));
}


} // namespace CGALi


template <class K>
inline
Object
intersection(
    const Iso_rectangle_2<K> &irect1,
    const Iso_rectangle_2<K> &irect2)
{
  return typename K::Intersect_2()(irect1, irect2);
}


template <class K>
inline bool
do_intersect(
    const Iso_rectangle_2<K> &irect1,
    const Iso_rectangle_2<K> &irect2)
{
    Object obj(intersection(irect1, irect2));
    typename K::Iso_rectangle_2 irect;
    return (assign(irect, obj));
}
CGAL_END_NAMESPACE

#endif
