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
// file          : src/Bbox_2_intersections.C
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#include <CGAL/Cartesian.h>
#include <CGAL/Bbox_2_Line_2_intersection.h>
#include <CGAL/Ray_2_Bbox_2_intersection.h>


CGAL_BEGIN_NAMESPACE

template <>
Bbox_2_Line_2_pair< Cartesian<double> >::Bbox_2_Line_2_pair(
    Bbox_2 const *bbox, Line_2< Cartesian<double> > const *line)
{
    _bbox = bbox;
    _line = *line;
    _known = false;
}

CGAL_END_NAMESPACE




CGAL_BEGIN_NAMESPACE

template <>
Ray_2_Bbox_2_pair< Cartesian<double> >::
Ray_2_Bbox_2_pair(
            Ray_2< Cartesian<double> > const *ray,
            Bbox_2 const *box)
{
    _known = false;
    _box = box;
    _ref_point = ray->start();
    _dir = ray->direction().to_vector();
    _min = 0;
}

CGAL_END_NAMESPACE


