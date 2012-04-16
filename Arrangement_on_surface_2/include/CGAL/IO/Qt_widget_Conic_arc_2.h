// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Ron Wein  <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_QT_WIDGET_CONIC_ARC_2_H
#define CGAL_QT_WIDGET_CONIC_ARC_2_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <list>

namespace CGAL {

/*!
 * Draw an x-monotone conic arc.
 */
template <class ConicArc>
Qt_widget& operator<< (Qt_widget& ws,
                       const _Conic_x_monotone_arc_2<ConicArc>& cv)
{
  // Get the co-ordinates of the curve's source and target.
  const double  sx = CGAL::to_double(cv.source().x());
  const double  sy = CGAL::to_double(cv.source().y());
  const double  tx = CGAL::to_double(cv.target().x());
  const double  ty = CGAL::to_double(cv.target().y());

  if (cv.orientation() == COLLINEAR)
  {
    // The curve is a segment - simply draw it.
    ws.get_painter().drawLine(ws.x_pixel(sx), ws.y_pixel(sy),
                              ws.x_pixel(tx), ws.y_pixel(ty));
    return (ws); 
  }

  // Draw a curves conic arc: As the arc is x-monotone, its source and its 
  // target has the extreme x-coordinates.
  const bool    is_source_left = (sx < tx);
  const int     x_min = is_source_left ? ws.x_pixel(sx) : ws.x_pixel(tx);
  const int     x_max = is_source_left ? ws.x_pixel(tx) : ws.x_pixel(sx);
  const int     n = x_max - x_min + 1;

  if (n <= 0)
    return (ws);
  
  typedef std::pair<double, double>    App_point_2;
  int                                  i;
  
  App_point_2  *pts = new App_point_2 [n + 1];
  cv.polyline_approximation (n, pts);

  ws.get_painter().moveTo (ws.x_pixel(pts[0].first),
                           ws.y_pixel(pts[0].second));
  for (i = 1; i <= n; i++)
  {
    ws.get_painter().lineTo (ws.x_pixel(pts[i].first),
                             ws.y_pixel(pts[i].second));
  }
  delete[] pts;
  
  return (ws);
}

/*!
 * Draw a conic arc.
 */
template <class Rat_kernel, class Alg_kernel, class Nt_traits>
Qt_widget& operator<< 
  (Qt_widget& ws, 
   const typename Arr_conic_traits_2<Rat_kernel,
                                     Alg_kernel,
                                     Nt_traits>::Curve_2& cv)
{
  typedef Arr_conic_traits_2<Rat_kernel,
                             Alg_kernel,
                             Nt_traits>                Conic_traits_2;
  typedef typename Conic_traits_2::X_monotone_curve_2  X_monotone_conic_arc_2;


  // Break the arc into x-monotone sub-curves and draw each one separately.
  Conic_traits_2                                              traits;
  std::list<X_monotone_conic_arc_2>                           x_arcs;
  typename std::list<X_monotone_conic_arc_2>::const_iterator  x_iter;

  traits.curve_make_x_monotone (cv,
				std::back_inserter (x_arcs));

  for (x_iter = x_arcs.begin(); x_iter != x_arcs.end(); ++x_iter)
    ws << *x_iter;

  return (ws); 
}

} //namespace CGAL

#endif
