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
// maintainer    : Efi Fogel <efif@post.tau.ac.il>
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <danha@post.tau.ac.il>)
//
// ======================================================================

#ifdef CGAL_ARR_POLYLINE_TRAITS_H
#ifndef CGAL_ARR_POLYLINE_TRAITS_WINDOW_STREAM_H   
#define CGAL_ARR_POLYLINE_TRAITS_WINDOW_STREAM_H  

#include <CGAL/Segment_2.h>
#include <CGAL/IO/Window_stream.h>    

CGAL_BEGIN_NAMESPACE

template <class Segment_traits_>
Window_stream& operator<< (Window_stream& ws,
			   const Polyline_2<Segment_traits_>& pl)
{ 
  typedef Polyline_2<Segment_traits_>          Curve_2;
  typedef typename Curve_2::const_iterator     Points_iterator;
  typedef typename Curve_2::Segment_2          Segment_2;

  Points_iterator   its = pl.begin();

  // Disregard empty polylines:
  if (its == pl.end())
    return (ws);

  // Draw the first point.
  ws << (*its);

  // Draw each segment of the polyline.
  Points_iterator   itt = pl.begin();
  itt++;

  while (itt != pl.end())
  {  
    ws << Segment_2(*its, *itt);
    its++; itt++;
  }

  // Now (*its) is the last polyline point -- draw it as well.
  ws << (*its);

  return (ws);
}

CGAL_END_NAMESPACE

#endif
#endif 









