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
// release       : $CGAL_Revision: CGAL-2.3-I-83 $
// release_date  : $CGAL_Date: 2001/07/17 $
//
// file          : include/CGAL/sweep_to_produce_planar_map_subcurves.h
// package       : Arrangement (2.10)
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

#ifndef SWEEP_TO_PRODUCE_PLANAR_MAP_POINTS_H
#define SWEEP_TO_PRODUCE_PLANAR_MAP_POINTS_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_SWEEP_CURVES_TO_SUBCURVES_H
#include <CGAL/Sweep_curves_to_subcurves.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Curve_iterator, class Traits, class Container>
void  sweep_to_produce_planar_map_points(Curve_iterator curves_begin, 
                                         Curve_iterator curves_end,  
                                         Traits& traits, 
                                         Container &points,
                                         bool overlapping = false)
{
  Sweep_curves_to_subcurves<Curve_iterator, Traits>  sweep_line;
  
  sweep_line.sweep_curves_to_points(curves_begin, curves_end, 
                                    points, overlapping);
}

/*template <class Curve_iterator, class PM, class Notifier>
  void  sweep_curves_to_pm(Curve_iterator curves_begin, 
  Curve_iterator curves_end, PM &result, Notifier* notifier)
  {
  Sweep_line<Curve_iterator, PM, Notifier>  sweep_line;
  
  sweep_line.sweep_curves_to_pm(curves_begin, curves_end, result, notifier);
  }*/

CGAL_END_NAMESPACE

#endif








