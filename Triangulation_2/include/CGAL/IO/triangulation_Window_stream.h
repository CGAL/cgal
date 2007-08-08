// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Olivier Devillers
//                 Andreas Fabri
//                 Monique Teillaud
//                 Mariette Yvinec


#ifdef CGAL_TRIANGULATION_2_H
#ifndef CGAL_WINDOW_STREAM_TRIANGULATION_2_H
#define CGAL_WINDOW_STREAM_TRIANGULATION_2_H

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds>
Window_stream&
operator<<(Window_stream& os,  const Triangulation_2<Gt, Tds> &t)
{
  return t.draw_triangulation(os);
}
CGAL_END_NAMESPACE
#endif // CGAL_WINDOW_STREAM_TRIANGULATION_2_H
#endif // CGAL_TRIANGULATION_2_H

#ifdef CGAL_DELAUNAY_TRIANGULATION_2_H
#ifndef CGAL_WINDOW_STREAM_DELAUNAY_TRIANGULATION_2_H
#define CGAL_WINDOW_STREAM_DELAUNAY_TRIANGULATION_2_H
CGAL_BEGIN_NAMESPACE
template < class Gt, class Tds >
Window_stream&
operator<<(Window_stream& os,  const Delaunay_triangulation_2<Gt,Tds> &dt)
{
   return dt.draw_triangulation(os);
}
CGAL_END_NAMESPACE
#endif // CGAL_WINDOW_STREAM_DELAUNAY_TRIANGULATION_2_H
#endif // CGAL_DELAUNAY_TRIANGULATION_2_H

#ifdef CGAL_CONSTRAINED_TRIANGULATION_2_H
#ifndef CGAL_WINDOW_STREAM_CONSTRAINED_TRIANGULATION_2_H
#define CGAL_WINDOW_STREAM_CONSTRAINED_TRIANGULATION_2_H
CGAL_BEGIN_NAMESPACE
template < class Gt, class Tds>
Window_stream&
operator<<(Window_stream& os,  const Constrained_triangulation_2<Gt,Tds> &t)
{
  return t.draw_triangulation(os);
}
CGAL_END_NAMESPACE
#endif // CGAL_WINDOW_STREAM_CONSTRAINED_TRIANGULATION_2_H
#endif // CGAL_CONSTRAINED_TRIANGULATION_2_H


#ifdef CGAL_WEIGHTED_POINT_H
#ifndef CGAL_WINDOW_STREAM_WEIGHTED_POINT_H
#define CGAL_WINDOW_STREAM_WEIGHTED_POINT_H
CGAL_BEGIN_NAMESPACE
template < class Point, class We >
Window_stream&
operator<<(Window_stream& os, const Weighted_point< Point, We > &p)
{
  double cx = to_double(p.point().x()),
         cy = to_double(p.point().y()),
         rr = to_double(p.weight());

  os<<p.point();
  os.draw_circle(cx, cy , std::sqrt(rr));
  return os;
}

template < class Point, class We >
Window_stream& operator>>(Window_stream &os, Weighted_point< Point, We > &wp)
{
  double cx, cy, x1, y1;
  os.read_mouse(cx,cy);
  os.read_mouse_circle(cx,cy, x1, y1);
  //os.read_mouse(x1,y1);
  Point center(cx,cy);

  We sr = We(std::sqrt( CGAL_NTS square(cx-x1)+ 
				 CGAL_NTS square(cy-y1) ) );

  os.draw_circle(cx, cy , sr);
  wp = Weighted_point< Point, We >(center, CGAL_NTS square(sr));
  return os;
}

CGAL_END_NAMESPACE
#endif // CGAL_WINDOW_STREAM_WEIGHTED_POINT_H
#endif // CGAL_WEIGHTED_POINT_H


#ifdef CGAL_REGULAR_TRIANGULATION_2_H
#ifndef CGAL_WINDOW_STREAM_REGULAR_TRIANGULATION_2_H
#define CGAL_WINDOW_STREAM_REGULAR_TRIANGULATION_2_H
CGAL_BEGIN_NAMESPACE
template < class Gt, class Tds >
Window_stream&
operator<<(Window_stream& os, Regular_triangulation_2<Gt,Tds> &t)
{
    return t.draw_triangulation(os);
}
CGAL_END_NAMESPACE
#endif // CGAL_WINDOW_STREAM_REGULAR_TRIANGULATION_2_H
#endif // CGAL_REGULAR_TRIANGULATION_2_H
