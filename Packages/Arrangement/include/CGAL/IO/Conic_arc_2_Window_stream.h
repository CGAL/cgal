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
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/Conic_arc_2_Window_stream.h
// package       : Arrangement (2.19)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// author(s)     : Ron Wein <wein@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CONIC_ARC_2_WINDOW_STREAM_H
#define CONIC_ARC_2_WINDOW_STREAM_H

#include <CGAL/IO/Window_stream.h>
#include <CGAL/Arrangement_2/Conic_arc_2.h>

CGAL_BEGIN_NAMESPACE

template <class Kernel>
Window_stream& operator<<(Window_stream& ws,
                          const Conic_arc_2<Kernel>& cv)
{
  // Get the co-ordinates of the curve's source and target.
  double sx = CGAL::to_double(cv.source().x()),
         sy = CGAL::to_double(cv.source().y()),
         tx = CGAL::to_double(cv.target().x()),
         ty = CGAL::to_double(cv.target().y());

  if (cv.is_segment())
  {
    // The curve is a segment - simply draw it.
    ws << leda_segment(sx, sy, tx, ty);
  } 
  // The arc is circular
  else
  {
    // If the curve is monotone, than its source and its target has the
    // extreme x co-ordinates on this curve.
    if (cv.is_x_monotone())
    {
      bool     is_source_left = (sx < tx);
      int      x_min = is_source_left ? ws.real_to_pix(sx) : 
	                                ws.real_to_pix(tx);
      int      x_max = is_source_left ? ws.real_to_pix(tx) :
                                        ws.real_to_pix(sx);
      double   prev_x = ws.pix_to_real(x_min);
      double   prev_y = is_source_left ? sy : ty;
      double   end_x = is_source_left ? tx : sx;
      double   end_y = is_source_left ? ty : sy;
      double   curr_x, curr_y;
      int      x;

      typename Conic_arc_2<Kernel>::Point_2 ps[2];
      int                      nps;

      for (x = x_min + 1; x < x_max; x++)
      {
	curr_x = ws.pix_to_real(x);
	nps = cv.get_points_at_x
          (Conic_arc_2<Kernel>::Point_2(typename Kernel::FT(curr_x), 0), ps);
	if (nps == 1)
	{
	  curr_y = CGAL::to_double(ps[0].y());
	  ws << leda_segment(prev_x, prev_y, curr_x, curr_y);
	  prev_x = curr_x;
	  prev_y = curr_y;
	  
	}
      }

      ws << leda_segment(prev_x, prev_y, end_x, end_y);
    }
    else
    {
      // We should never reach here.
      CGAL_assertion(false);
    }
  }

  return (ws); 
}

CGAL_END_NAMESPACE

#endif
