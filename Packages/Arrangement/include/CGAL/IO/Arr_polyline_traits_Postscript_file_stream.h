// Copyright (c) 2001  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>

#ifdef CGAL_ARR_POLYLINE_TRAITS_H
#ifndef CGAL_ARR_POLYLINE_TRAITS_POSTSCRIPT_FILE_STREAM_H   
#define CGAL_ARR_POLYLINE_TRAITS_POSTSCRIPT_FILE_STREAM_H  

#include <CGAL/Segment_2.h> 
#include <CGAL/IO/Postscript_file_stream.h>

CGAL_BEGIN_NAMESPACE

template <class Segment_traits_>
Postscript_file_stream& operator<< (Postscript_file_stream& ps,
				    const Polyline_2<Segment_traits_>& pl)
{
  typedef Polyline_2<Segment_traits_>          Curve_2;
  typedef typename Curve_2::const_iterator     Points_iterator;
  typedef typename Curve_2::Segment_2          Segment_2;

  Points_iterator   its = pl.begin();

  // Disregard empty polylines:
  if (its == pl.end())
    return (ps);

  // Draw the first point.
  ps << (*its);

  // Draw each segment of the polyline.
  Points_iterator   itt = pl.begin();
  itt++;

  while (itt != pl.end())
  {  
    ps << Segment_2(*its, *itt);
    its++; itt++;
  }

  // Now (*its) is the last polyline point -- draw it as well.
  ps << (*its);

  return (ps);
}

CGAL_END_NAMESPACE

#endif
#endif









