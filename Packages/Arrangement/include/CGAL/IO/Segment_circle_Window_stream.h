#ifndef SEGMENT_CIRCLE_WINDOW_STREAM_H
#define SEGMENT_CIRCLE_WINDOW_STREAM_H

#include <CGAL/IO/Window_stream.h>

#include <CGAL/Segment_circle_2.h>

CGAL_BEGIN_NAMESPACE

template <class NT>
Window_stream& operator<<(Window_stream& os,
                          const Segment_circle_2<NT>& cv)
{
  double sx = CGAL::to_double(cv.source().x()),
         sy = CGAL::to_double(cv.source().y()),
         tx = CGAL::to_double(cv.target().x()),
         ty = CGAL::to_double(cv.target().y());

  if (cv.is_segment())
  {
    W << leda_segment(sx, sy, tx, ty);
  } 
  // The arc is circular
  else
  {
    // We need a middle point on the curve for the leda draw_arc
    // function which draws a circular arc given three points on it.
    
    double px,py; // middle point coordinates
    if (cv.is_x_monotone())
    {
      // an x-monotone circular arc 
      // the middle point is the one with average x value of endpoints.
      Segment_circle_2<NT>::Point ps[2];
      NT middle_x = (sx+tx)/2;
      cv.get_points_at_x(middle_x, ps);
      px = CGAL::to_double(middle_x);
      py = CGAL::to_double(ps[0].y());
     }
    else
    {
      // a non x-monotone circular arc
      // we use the rightmost or leftmost point as a middle point
      Segment_circle_2<NT>::Point ps[2];
      // there might be two tangency points but we care for one only
      cv.horizontal_tangency_points (ps);
      px = CGAL::to_double(ps[0].x());
      py = CGAL::to_double(ps[0].y());
    }
    
    os.draw_arc(leda_point(sx,sy), leda_point(px,py), leda_point(tx,ty));
  }

  return os; 
}

CGAL_END_NAMESPACE

#endif
