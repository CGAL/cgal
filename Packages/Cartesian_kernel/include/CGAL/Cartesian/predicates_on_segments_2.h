// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/predicates_on_segments_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Susan Hert
// coordinator   : INRIA Sophia-Antipolis
//
// ======================================================================

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
