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
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/sweep_to_produce_subcurves_2.h
// package       : arr (1.87)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef SWEEP_TO_PRODUCE_PLANAR_MAP_SUBCURVES_H
#define SWEEP_TO_PRODUCE_PLANAR_MAP_SUBCURVES_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_SWEEP_CURVES_TO_SUBCURVES_H
#include <CGAL/Sweep_line_2/Sweep_curves_to_subcurves_2.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Curve_iterator, class SweepLineTraits_2, class OutpoutIterator>
void  sweep_to_produce_subcurves_2(Curve_iterator curves_begin, 
                                   Curve_iterator curves_end,  
                                   SweepLineTraits_2& traits, 
                                   OutpoutIterator subcurves,
                                   bool overlapping = false)
{
  Sweep_curves_to_subcurves_2<Curve_iterator, SweepLineTraits_2>  
                                                        sweep_line(&traits);
  
  sweep_line.sweep_curves_to_subcurves(curves_begin, curves_end, 
                                       subcurves, overlapping);
}

CGAL_END_NAMESPACE

#endif




