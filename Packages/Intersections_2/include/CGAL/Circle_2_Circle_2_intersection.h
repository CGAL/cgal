
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
// file          : include/CGAL/Circle_2_Circle_2_intersection.h
// source        : intersection_2_3.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : MPI, Saarbruecken
//
// ============================================================================


#ifndef CGAL_CIRCLE_2_CIRCLE_2_INTERSECTION_H
#define CGAL_CIRCLE_2_CIRCLE_2_INTERSECTION_H

#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
bool
do_intersect(const Circle_2<R> &tr1, const Circle_2<R>&tr2);

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Circle_2_Circle_2_intersection.C>
#endif
#endif
