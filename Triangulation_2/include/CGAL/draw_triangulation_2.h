// Copyright(c) 2018  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_DRAW_T2_H
#define CGAL_DRAW_T2_H

#include <CGAL/Drawing_functor.h>
#include <CGAL/Graphic_buffer.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/license/Triangulation_2.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_2.h>

namespace CGAL {

namespace draw_function_for_t2 {

template <typename BufferType = float, class T2, class DrawingFunctor>
void compute_face(typename T2::Finite_faces_iterator fh,
                  const DrawingFunctor &drawing_functor,
                  CGAL::GraphicBuffer<BufferType> &graphic_buffer,
                  const T2 *t2) {

  if (!drawing_functor.draw_face(*t2, fh)) {
    return;
  }

  if (drawing_functor.colored_face(*t2, fh) && drawing_functor.face_color) {
    CGAL::IO::Color c = drawing_functor.face_color(*t2, fh);
    graphic_buffer.face_begin(c);

  } else {
    graphic_buffer.face_begin();
  }

  graphic_buffer.add_point_in_face(fh->vertex(0)->point());
  graphic_buffer.add_point_in_face(fh->vertex(1)->point());
  graphic_buffer.add_point_in_face(fh->vertex(2)->point());

  graphic_buffer.face_end();
}

template <typename BufferType = float, class T2, class DrawingFunctor>
void compute_edge(typename T2::Finite_edges_iterator eh,
                  const DrawingFunctor &drawing_functor, const T2 *t2,
                  CGAL::GraphicBuffer<BufferType> &graphic_buffer) {

  if (!drawing_functor.draw_edge(*t2, eh)) {
    return;
  }

  if (drawing_functor.colored_edge(*t2, eh) && drawing_functor.edge_color) {

    CGAL::IO::Color c = drawing_functor.edge_color(*t2, eh);
    graphic_buffer.add_segment(
        eh->first->vertex(eh->first->cw(eh->second))->point(),
        eh->first->vertex(eh->first->ccw(eh->second))->point(), c);
  }

  else {
    graphic_buffer.add_segment(
        eh->first->vertex(eh->first->cw(eh->second))->point(),
        eh->first->vertex(eh->first->ccw(eh->second))->point());
  }
}

template <typename BufferType = float, class T2, class DrawingFunctor>
void compute_vertex(typename T2::Vertex_handle vh,
                    CGAL::GraphicBuffer<BufferType> &graphic_buffer,
                    const DrawingFunctor &drawing_functor, const T2 *t2) {

  if (!drawing_functor.draw_vertex(*t2, vh)) {
    return;
  }

  if (drawing_functor.colored_vertex(*t2, vh) && drawing_functor.vertex_color) {
    CGAL::IO::Color c = drawing_functor.vertex_color(*t2, vh);
    graphic_buffer.add_point(vh->point(), c);
  } else {
    graphic_buffer.add_point(vh->point());
  }
}

template <typename BufferType = float, class T2, class DrawingFunctor>
void compute_elements(CGAL::GraphicBuffer<BufferType> &graphic_buffer,
                      const T2 *t2, const DrawingFunctor &drawing_functor,
                      bool m_nofaces = false) {

  if (!m_nofaces) {
    for (typename T2::Finite_faces_iterator it = t2->finite_faces_begin();
         it != t2->finite_faces_end(); ++it) {
      compute_face(it, drawing_functor, graphic_buffer, t2);
    }
  }

  for (typename T2::Finite_edges_iterator it = t2->finite_edges_begin();
       it != t2->finite_edges_end(); ++it) {
    compute_edge(it, drawing_functor, t2, graphic_buffer);
  }

  for (typename T2::Finite_vertices_iterator it = t2->finite_vertices_begin();
       it != t2->finite_vertices_end(); ++it) {
    compute_vertex(it, graphic_buffer, drawing_functor, t2);
  }
}

} // namespace draw_function_for_t2

template <typename BufferType = float, class T2, class DrawingFunctor>
void add_in_graphic_buffer_t2(CGAL::GraphicBuffer<BufferType> &graphic_buffer,
                              const DrawingFunctor &drawing_functor,
                              const T2 *at2 = nullptr, bool m_nofaces = false) {
  if (at2 != nullptr) {
    draw_function_for_t2::compute_elements(graphic_buffer, at2, drawing_functor,
                                           m_nofaces);
  }
}

// Specialization of draw function.
#define CGAL_T2_TYPE CGAL::Triangulation_2<Gt, Tds>

template <class Gt, class Tds, class DrawingFunctor>
void draw(const CGAL_T2_TYPE &at2, const DrawingFunctor &drawingfunctor,
          const char *title = "Triangulation_2 Basic Viewer",
          bool nofill = false) {

  CGAL::GraphicBuffer<float> buffer;
  add_in_graphic_buffer_t2(buffer, drawingfunctor, &at2, nofill);
  draw_buffer(buffer);
}

template <class Gt, class Tds>
void draw(const CGAL_T2_TYPE &at2,
          const char *title = "Triangulation_2 Basic Viewer",
          bool nofill = false) {

  CGAL::GraphicBuffer<float> buffer;
  Drawing_functor<CGAL_T2_TYPE, typename CGAL_T2_TYPE::Vertex_handle,
                  typename CGAL_T2_TYPE::Finite_edges_iterator,
                  typename CGAL_T2_TYPE::Finite_faces_iterator>
      drawingFunctor;

  drawingFunctor.colored_face =
      [](const CGAL_T2_TYPE &at3,
         const typename CGAL_T2_TYPE::Finite_faces_iterator fh) -> bool {
    return true;
  };

  drawingFunctor.face_color =
      [](const CGAL_T2_TYPE &at2,
         const typename CGAL_T2_TYPE::Finite_faces_iterator fh)
      -> CGAL::IO::Color {
    CGAL::Random random((unsigned int)(std::size_t)(&*fh));
    return get_random_color(random);
  };

  draw(at2, drawingFunctor, title, nofill);
}

#undef CGAL_T2_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_T2_H
