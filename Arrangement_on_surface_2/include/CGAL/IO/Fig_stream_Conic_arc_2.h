// Copyright (c) 2006,2007,2009,2011 Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
// Author(s)     : Ron Wein           <wein@post.tau.ac.il>

#ifndef CGAL_FIG_STREAM_CONIC_ARC_2_H
#define CGAL_FIG_STREAM_CONIC_ARC_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/IO/Fig_stream.h>
#include <list>

/*!
 * Write an x-monotone conic arc to a FIG stream.
 */
template <class Conic_traits>
static void _write_x_monotone_conic_arc 
   (CGAL::Fig_stream<typename Conic_traits::Alg_kernel>& fs,
    const typename Conic_traits::X_monotone_curve_2&     cv)
{
  typedef typename Conic_traits::Alg_kernel  Alg_kernel;
  typedef typename Conic_traits::Algebraic   Algebraic;
  typedef typename Conic_traits::Point_2     Alg_point_2;
  typedef typename Alg_kernel::Segment_2     Alg_segment_2;
  
  if (cv.orientation() == CGAL::COLLINEAR)
  {
    // In case of a linear segment:
    Alg_segment_2   seg = Alg_segment_2 (cv.source(), cv.target());
    fs << seg;
  } 
  else if (CGAL::compare (cv.r(), cv.s()) == CGAL::EQUAL &&
	   CGAL::sign (cv.t()) == CGAL::ZERO)
  {
    // In case of a circular arc:
    Algebraic    x_mid = (cv.source().x() + cv.target().x()) / 2;
    Alg_point_2  q = Alg_point_2(x_mid, 0);
    Alg_point_2  p = cv.point_at_x (q);

    fs.write_circular_arc (cv.source(), p, cv.target());
  }
  else
  {
    // Represent the arc as a spline with 5 control points.
    Algebraic   x;
    Alg_point_2 q;
    Alg_point_2 cps[5];
    int         i;

    cps[0] = cv.source();
    for (i = 1; i <= 3; i++)
    {
      x = (cv.source().x()*(4 - i) + cv.target().x()*i) / 4;
      
      q = Alg_point_2(x, 0);
      cps[i] = cv.point_at_x (q);      
    }
    cps[4] = cv.target();

    fs.write_spline (cps + 0, cps + 5, 1.0);
  }

  return;
}

/*!
 * Write a conic arc to a FIG stream.
 */
template <class Conic_traits>
void write_conic_arc
   (CGAL::Fig_stream<typename Conic_traits::Alg_kernel>& fs,
    const typename Conic_traits::Curve_2&                cv)
{
  typedef typename Conic_traits::X_monotone_curve_2   Conic_arc_2;

  // Subdivide the arc into x-monotone sub-arcs.
  Conic_traits           traits;
  std::list<Conic_arc_2> xcvs;

  traits.make_x_monotone_2_object() (cv,
				     std::back_inserter(xcvs));

  // Write the x-monotone sub-arcs.
  typename std::list<Conic_arc_2>::iterator xit;
    
  for (xit = xcvs.begin(); xit != xcvs.end(); ++xit)
    _write_x_monotone_conic_arc<Conic_traits> (fs, *xit);

  return;
}

#endif
