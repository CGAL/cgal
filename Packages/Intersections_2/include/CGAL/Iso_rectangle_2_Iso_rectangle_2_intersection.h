
// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision:  $
// release_date  : $CGAL_Date:  $
//
// file          : include/CGAL/Iso_rectangle_2_Iso_rectangle_2_intersection.h
// source        : intersection_2_2.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


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
        return K::Object_2(); 
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
  return CGALi::intersection(irect1, irect2, K());
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
