// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-62 $
// release_date  : $CGAL_Date: 2001/05/11 $
//
// file          : include/CGAL/IO/Arr_polyline_traits_Window_stream.h
// package       : Arrangement (1.82)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================

#ifdef CGAL_ARR_POLYLINE_TRAITS_H
#ifndef CGAL_ARR_POLYLINE_TRAITS_WINDOW_STREAM_H   
#define CGAL_ARR_POLYLINE_TRAITS_WINDOW_STREAM_H  

#include <CGAL/Segment_2.h>
#include <CGAL/IO/Window_stream.h>    

CGAL_BEGIN_NAMESPACE

template <class Curve>
Window_stream& operator<<(Window_stream& W, const Curve& cv)        
{ 
  typedef typename Curve::value_type           Point;
  typedef typename Curve::iterator             Curve_iter;
  typedef typename Point::R                    R;
  typedef CGAL::Segment_2<R>                   Segment;
  typedef typename Curve::const_iterator       Points_iterator;

  W << *(cv.begin());
  
  Points_iterator points_iter;
  for (points_iter = cv.begin(); points_iter != cv.end(); ) {
    //for (unsigned int i = 0; i < cv.size() - 1; i++){
    Points_iterator source_point = points_iter;
    
    Points_iterator target_point =  (++points_iter);
    
    if (target_point == cv.end())
      break;
    
    W << Segment(*source_point, *target_point);
  }

  W << *(--points_iter);

  return W;
}

CGAL_END_NAMESPACE

#endif
#endif 









