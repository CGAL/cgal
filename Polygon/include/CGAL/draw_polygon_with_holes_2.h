// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_DRAW_POLYGON_WITH_HOLES_2_H
#define CGAL_DRAW_POLYGON_WITH_HOLES_2_H

#include <CGAL/Drawing_functor.h>
#include <CGAL/Graphic_buffer.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
\ingroup PkgDrawPolygonWithHoles2

opens a new window and draws `aph`, an instance of the `CGAL::Polygon_with_holes_2` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam PH an instance of the `CGAL::Polygon_with_holes_2` class.
\param aph the polygon with holes to draw.

*/
template<class PH>
void draw(const PH& aph);

} /* namespace CGAL */
#endif

#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Random.h>

namespace CGAL
{

namespace draw_function_for_ph2_with_holes {

template <typename BufferType = float, class P2>
void compute_one_loop_elements(const typename P2::General_polygon_2& p, Graphic_buffer<BufferType> &graphic_buffer, bool hole)
{
  
  if (hole)
  { graphic_buffer.add_point_in_face(p.vertex(p.size()-1)); }

  typename P2::General_polygon_2::Vertex_const_iterator prev;
  for (typename P2::General_polygon_2::Vertex_const_iterator i=p.vertices_begin();
        i!=p.vertices_end(); ++i)
  {
    graphic_buffer.add_point(*i);         // Add vertex
    if (i!=p.vertices_begin())
    { graphic_buffer.add_segment(*prev, *i); } // Add segment with previous point
    graphic_buffer.add_point_in_face(*i); // Add point in face
    prev=i;
  }

  // Add the last segment between the last point and the first one
  graphic_buffer.add_segment(*prev, *(p.vertices_begin()));
}

template <typename BufferType = float, class P2, class DrawingFunctor>
void compute_elements(const P2& p2, Graphic_buffer<BufferType> &graphic_buffer, 
                      const DrawingFunctor &m_drawing_functor)
{

  if (p2.outer_boundary().is_empty()) return;

  // TODO: use face_color after adding a handler if exists.
  CGAL::IO::Color c(75,160,255);
  graphic_buffer.face_begin(c);

  compute_one_loop_elements<BufferType, P2>(p2.outer_boundary(), graphic_buffer, false);

  for (typename P2::Hole_const_iterator it=p2.holes_begin(); it!=p2.holes_end(); ++it)
  {
    compute_one_loop_elements<BufferType, P2>(*it, graphic_buffer, true);
    graphic_buffer.add_point_in_face(p2.outer_boundary().vertex(p2.outer_boundary().size()-1));
  }

  graphic_buffer.face_end();
}

} // draw_function_for_ph2


template <typename BufferType = float, class P2, class DrawingFunctor>
void add_in_graphic_buffer(const P2 &p2, CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                             const DrawingFunctor &m_drawing_functor ) {
  draw_function_for_ph2_with_holes::compute_elements(p2, graphic_buffer, m_drawing_functor);
}

template <typename BufferType = float, class P2>
void add_in_graphic_buffer(const P2 &p2, CGAL::Graphic_buffer<BufferType> &graphic_buffer) {

  // TODO: use colord_face and face_color if a handler exits.
  Drawing_functor<P2, typename P2::Hole_const_iterator,
                  typename P2::Hole_const_iterator,
                  typename P2::Hole_const_iterator>
      drawing_functor;

  add_in_graphic_buffer(p2, graphic_buffer, drawing_functor);
}

// Specialization of draw function.
#define CGAL_P2_WITH_HOLES_TYPE CGAL::Polygon_with_holes_2<T, C>


template<class T, class C, typename BufferType = float, class DrawingFunctor>
void draw(const CGAL_P2_WITH_HOLES_TYPE& ap2, const DrawingFunctor &drawing_functor)
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(ap2, buffer, drawing_functor);
  draw_buffer(buffer);
}

template<class T, class C, typename BufferType = float>
void draw(const CGAL_P2_WITH_HOLES_TYPE& ap2,
          const char* title="Polygon_with_holes_2 Basic Viewer")
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(ap2, buffer);
  draw_buffer(buffer);
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYGON_WITH_HOLES_2_H
