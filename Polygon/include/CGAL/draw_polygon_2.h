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

#ifndef CGAL_DRAW_POLYGON_2_H
#define CGAL_DRAW_POLYGON_2_H

#include <CGAL/Drawing_functor.h>
#include <CGAL/Graphic_buffer.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
\ingroup PkgDrawPolygon2

opens a new window and draws `ap`, an instance of the `CGAL::Polygon_2` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam P an instance of the `CGAL::Polygon_2` class.
\param ap the polygon to draw.

*/
template<class P>
void draw(const P& ap);

} /* namespace CGAL */

#endif

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Polygon_2.h>
#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Random.h>

namespace CGAL {

namespace draw_function_for_p2 {

template <typename BufferType = float, class P2>
void compute_elements(CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                      const P2 *p2) {

  typedef typename P2::Point_2 Point;

  if (p2->is_empty())
    return;

  Point prev = p2->vertex(p2->size() - 1);

  CGAL::IO::Color c(75, 160, 255);
  graphic_buffer.face_begin(c);

  for (typename P2::Vertex_const_iterator i = p2->vertices_begin();
       i != p2->vertices_end(); ++i) {
    graphic_buffer.add_point(*i);         // Add vertex
    graphic_buffer.add_segment(prev, *i); // Add segment with previous point
    graphic_buffer.add_point_in_face(*i); // Add point in face
    prev = *i;
  }

  graphic_buffer.face_end();
}

} // namespace draw_function_for_p2

template <typename BufferType = float, class P2>
void add_in_graphic_buffer(CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                              const P2 *p2 = nullptr) {
  if (p2 != nullptr) {
    draw_function_for_p2::compute_elements(graphic_buffer, p2);
  }
}

// Specialization of draw function.
#define CGAL_P2_TYPE CGAL::Polygon_2<T, C>

template <class T, class C>
void draw(const CGAL_P2_TYPE &ap2,
          const char *title = "Polygon_2 Basic Viewer") {

  // Drawing_functor<CGAL_P2_TYPE, typename CGAL_P2_TYPE::Vertex_const_handle,
  //                 typename CGAL_P2_TYPE::Halfedge_const_handle,
  //                 typename CGAL_P2_TYPE::Face_const_handle>
  //     drawingFunctor;

  CGAL::Graphic_buffer<float> buffer;
  add_in_graphic_buffer(buffer, &ap2);
  draw_buffer(buffer);
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYGON_2_H
