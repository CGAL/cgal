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

template < class R >
Comparison_result
compare_slopes(const SegmentC2<R CGAL_CTAG> &s1,
               const SegmentC2<R CGAL_CTAG> &s2)
{
   return compare_slopesC2(s1.source().x(), s1.source().y(),
                           s1.target().x(), s1.target().y(),
                           s2.source().x(), s2.source().y(),
}


CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PREDICATES_ON_SEGMENTS_2_H
