
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
// file          : include/CGAL/Circle_2_Circle_2_intersection.C
// source        : intersection_2_2.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================



#include <CGAL/Circle_2.h>
#include <CGAL/squared_distance_2_1.h>
CGAL_BEGIN_NAMESPACE

template <class R>
bool
do_intersect(const Circle_2<R> & circ1, const Circle_2<R>& circ2)
{
    typedef typename R::FT FT;
    FT sr1 = circ1.squared_radius();
    FT sr2 = circ2.squared_radius();
    FT squared_dist = squared_distance(circ1.center(), circ2.center());
    FT temp = sr1+sr2-squared_dist;
    return !(FT(4)*sr1*sr2 < temp*temp);
}

CGAL_END_NAMESPACE

