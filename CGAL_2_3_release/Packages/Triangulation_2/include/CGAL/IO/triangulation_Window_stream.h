// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/IO/triangulation_Window_stream.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Olivier Devillers
//                 Andreas Fabri
//                 Monique Teillaud
//                 Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


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
  os.draw_circle(cx, cy , CGAL_CLIB_STD::sqrt(rr));
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

  We sr = We(CGAL_CLIB_STD::sqrt( CGAL_NTS square(cx-x1)+ 
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
