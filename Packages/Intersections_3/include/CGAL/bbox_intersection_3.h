
// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : include/CGAL/bbox_intersection_3.h
// source        : web/intersection_3.fw
// author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_BBOX_INTERSECTION_3_H
#define CGAL_BBOX_INTERSECTION_3_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Object.h>


CGAL_BEGIN_NAMESPACE

extern Object
intersection_bl(const Bbox_3 &box,
        double lx1, double ly1, double lz1,
        double lx2, double ly2, double lz2,
        bool min_infinite, bool max_infinite);
CGAL_END_NAMESPACE


#endif // CGAL_BBOX_INTERSECTION_3_H
