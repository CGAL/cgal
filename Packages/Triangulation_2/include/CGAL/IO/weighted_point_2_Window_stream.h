#ifdef CGAL_WEIGHTED_POINT_2_H
#ifndef CGAL_WINDOW_STREAM_WEIGHTED_POINT_2_H
#define CGAL_WINDOW_STREAM_WEIGHTED_POINT_2_H

template < class Point, class We >
CGAL_Window_stream&
operator<<(CGAL_Window_stream& os,
           const CGAL_Weighted_point_2< Point, We > &p)
{
  double cx = CGAL_to_double(p.point().x()),
         cy = CGAL_to_double(p.point().y()),
         r = CGAL_to_double(p.weight());

  os<<p.point();
  //os.draw_circle(cx, cy , /*sqrt*/(r));
  return os;
}

template < class Point, class We >
CGAL_Window_stream& operator>>(CGAL_Window_stream &os,
                               CGAL_Weighted_point_2< Point, We > &wp)
{
  double cx, cy, x1, y1;
  os.read_mouse(cx,cy);
  os.read_mouse_circle(cx,cy, x1, y1);
  Point center(cx, cy);

  We sr = We (sqrt( (cx-x1)*(cx-x1)+(cy-y1)*(cy-y1) ) );

  os.draw_circle(cx, cy , sr );
  wp = CGAL_Weighted_point_2< Point, We >(center, sr);
  return os;
}

#endif // CGAL_WINDOW_STREAM_WEIGHTED_POINT_2_H
#endif // CGAL_WEIGHTED_POINT_2_H

