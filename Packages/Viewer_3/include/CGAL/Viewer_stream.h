// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Francois Rebufat <Francois.Rebufat@sophia.inria.fr>
#ifndef CGAL_VIEWER_STREAM_H
#define CGAL_VIEWER_STREAM_H

#include <CGAL/Point_3.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Viewer_3.h>

CGAL_BEGIN_NAMESPACE
// Transform a point_3 to a drawable_point_3 with Viewer attributs.
template <class R>
Drawable_point_3<Point_3<R> >*
convert_type(const Point_3<R> p, Viewer_3 &W)
{
  Drawable_point_3<Point_3<R> >* dp = new
    Drawable_point_3<Point_3<R> >(p,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(dp);
}

template <class R>
Drawable_point_3<std::vector<R> >*
convert_type(const std::vector<R> &p, Viewer_3 &W)
{
  Drawable_point_3<std::vector<R> >* dp = new
    Drawable_point_3<std::vector<R> >(p[0],p[1],p[2],W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(dp);
}



template<class R>
Drawable_segment_3<Segment_3<R> >*
convert_type(const Segment_3<R> &s, Viewer_3 &W)
{
   Drawable_segment_3<Segment_3<R> >* sp = new
     Drawable_segment_3<Segment_3<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}

template<class R>
Drawable_triangle_3<Triangle_3<R> >*
convert_type(const Triangle_3<R> &s, Viewer_3 &W)
{
   Drawable_triangle_3<Triangle_3<R> >* sp = new
     Drawable_triangle_3<Triangle_3<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}


template<class R>
Drawable_tetrahedron_3<Tetrahedron_3<R> >*
convert_type(const Tetrahedron_3<R> &s, Viewer_3 &W)
{
   Drawable_tetrahedron_3<Tetrahedron_3<R> >* sp = new
     Drawable_tetrahedron_3<Tetrahedron_3<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}

template<class R>
Drawable_line_3<Line_3<R> >*
convert_type(const Line_3<R> &s, Viewer_3 &W)
{
   Drawable_line_3<Line_3<R> >* sp = new
     Drawable_line_3<Line_3<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}


template<class R>
Drawable_ray_3<Ray_3<R> >*
convert_type(const Ray_3<R> &s, Viewer_3 &W)
{
   Drawable_ray_3<Ray_3<R> >* sp = new
     Drawable_ray_3<Ray_3<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}


template<class R>
Drawable_point_2<Point_2<R> >*
convert_type(const Point_2<R> &s, Viewer_3 &W)
{
   Drawable_point_2<Point_2<R> >* sp = new
     Drawable_point_2<Point_2<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}

template<class R>
Drawable_segment_2<Segment_2<R> >*
convert_type(const Segment_2<R> &s, Viewer_3 &W)
{
   Drawable_segment_2<Segment_2<R> >* sp = new
     Drawable_segment_2<Segment_2<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}

template<class R>
Drawable_line_2<Line_2<R> >*
convert_type(const Line_2<R> &s, Viewer_3 &W)
{
   Drawable_line_2<Line_2<R> >* sp = new
     Drawable_line_2<Line_2<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}

template<class R>
Drawable_ray_2<Ray_2<R> >*
convert_type(const Ray_2<R> &s, Viewer_3 &W)
{
   Drawable_ray_2<Ray_2<R> >* sp = new
     Drawable_ray_2<Ray_2<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}

template<class R>
Drawable_triangle_2<Triangle_2<R> >*
convert_type(const Triangle_2<R> &s, Viewer_3 &W)
{
   Drawable_triangle_2<Triangle_2<R> >* sp = new
     Drawable_triangle_2<Triangle_2<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}

template<class R>
Drawable_circle_2<Circle_2<R> >*
convert_type(const Circle_2<R> &s, Viewer_3 &W)
{
   Drawable_circle_2<Circle_2<R> >* sp = new
     Drawable_circle_2<Circle_2<R> >(s,W.get_color(1),W.get_style(),W.get_size(),W.get_precision());
  return(sp);
}

// ######## The Stream Manipulator ###############################
template<class Obj> class O_manip;
template< class Obj > Viewer_3& operator<<(Viewer_3& W, const
					   O_manip<Obj> &m);

template<class Obj> class O_manip {
#ifdef _MSC_VER
 public:
#endif
  Viewer_3& (*f)(Viewer_3&,Obj);
  Obj i;

public:
  O_manip(Viewer_3& (*ff)(Viewer_3&,Obj), Obj ii)
    : f(ff), i(ii){}

  friend Viewer_3& operator<< CGAL_NULL_TMPL_ARGS
  (Viewer_3& W, const O_manip<Obj>& m);
   
};


template<class Obj>
Viewer_3& operator<<(Viewer_3& W, const O_manip<Obj>& m)
{
  return m.f(W,m.i);
}

// ######## The Stream Manipulator END###############################

// functions used for precision
Viewer_3& precision(Viewer_3& W, Precision s)
{
  W.set_precision(s);
  return W;
}

O_manip<Precision> set_precision(Precision i)
{
  return O_manip<Precision>(&precision, i);
}

// functions used for colors
Viewer_3& color_1(Viewer_3& W, Color c)
{
  W.set_color(c,1);
  return W;
}

O_manip<Color> set_color_1(Color i)
{
  return O_manip<Color>(&color_1, i);
}
// --------------------
Viewer_3& color_2(Viewer_3& W, Color c)
{
  W.set_color(c,2);
  return W;
}

O_manip<Color> set_color_2(Color i)
{
  return O_manip<Color>(&color_2, i);
}
// functions used for size
Viewer_3& size_f(Viewer_3& W, Size c)
{
  W.set_size(c);
  return W;
}

O_manip<Size> set_size(Size i)
{
  return O_manip<Size>(&size_f, i);
}
// functions used for style
Viewer_3& style_f(Viewer_3& W, Style c)
{
  W.set_style(c);
  return W;
}

O_manip<Style> set_style(Style i)
{
  return O_manip<Style>(&style_f, i);
}

// }

// Stream for style and precision.
// Viewer_3&
// operator<<(Viewer_3& W, Style s)
// {
// W.set_style(s);
// return W;
// }

// Viewer_3&
// operator<<(Viewer_3& W, Precision s)
// {
// W.set_precision(s);
// return W;
// }
// Stream for a CGAL object that has a convert_type operation
template <class object>
Viewer_3&
operator<<(Viewer_3& W, const object &o)
{
W.add_drawable(convert_type(o,W));
return W;
}
CGAL_END_NAMESPACE

#endif //CGAL_VIEWER_STREAM_H

