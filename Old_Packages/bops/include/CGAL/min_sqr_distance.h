// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 01
//
// file          : include/CGAL/min_sqr_distance.h
// package       : bops (2.2)
// source        : include/CGAL/min_sqr_distance.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_MIN_SQR_DISTANCE_H
#define CGAL_MIN_SQR_DISTANCE_H

CGAL_BEGIN_NAMESPACE

template <class ForwardIterator, class I>
double minimal_square_distance2(ForwardIterator first,
			        ForwardIterator last,
			        I& T)
/*
   ----------------------------------------------------------
   Calculates the minimal square distance of a set of points
   in O(n^2) time and returns the (double) result;
   ----------------------------------------------------------
*/
;
 


template <class ForwardIterator, class I>
double minimal_square_distance(ForwardIterator first,
			       ForwardIterator last,
			       I& T)
/*
   ----------------------------------------------------------
   Calculates the minimal square distance of a set of points
   in O(n log(n)) time by doing a sweep-line algorithm and
   returns the (double) result;

   See:
     K.Hinrichs, J.Nievergelt, P.Schorn,
     "Plane-sweep solves the closest pair problem elegantly",
     Inform. Process. Lett., 26, 1988, p. 255-261.
   or
     Rolf Klein: "Algorithmische Geometrie", chapter 2.3.1
   ----------------------------------------------------------
*/
;

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/min_sqr_distance.C>
#endif

#endif // CGAL_MIN_SQR_DISTANCE_H
