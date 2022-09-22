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

#include <CGAL/Generic_functors.h>
#include <CGAL/GraphicBuffer.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/license/Triangulation_2.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_2.h>

namespace CGAL {

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorT2 {
  template <typename T2>
  static CGAL::IO::Color run(const T2 &,
                             const typename T2::Finite_faces_iterator fh) {
    CGAL::Random random((unsigned int)(std::size_t)(&*fh));
    return get_random_color(random);
  }
};

namespace draw_function_for_t2 {

template <typename BufferType = float, class T2, class ColorFunctor>
void compute_face(typename T2::Finite_faces_iterator fh,
                  GraphicBuffer<BufferType> &graphic_buffer, const T2 *t2,
                  const ColorFunctor &m_fcolor) {

  // TODO: change it to m_fcolor.face_color(*t2, fh)
  CGAL::IO::Color c = m_fcolor.run(*t2, fh);
  graphic_buffer.face_begin(c);

  graphic_buffer.add_point_in_face(fh->vertex(0)->point());
  graphic_buffer.add_point_in_face(fh->vertex(1)->point());
  graphic_buffer.add_point_in_face(fh->vertex(2)->point());

  graphic_buffer.face_end();
}

template <typename BufferType = float, class T2>
void compute_edge(typename T2::Finite_edges_iterator eh,
                  GraphicBuffer<BufferType> &graphic_buffer) {
  graphic_buffer.add_segment(
      eh->first->vertex(eh->first->cw(eh->second))->point(),
      eh->first->vertex(eh->first->ccw(eh->second))->point());
}

template <typename BufferType = float, class T2>
void compute_vertex(typename T2::Vertex_handle vh,
                    GraphicBuffer<BufferType> &graphic_buffer) {
  graphic_buffer.add_point(vh->point());
}

template <typename BufferType = float, class T2, class ColorFunctor>
void compute_elements(GraphicBuffer<BufferType> &graphic_buffer, const T2 *t2,
                      const ColorFunctor &m_color_functor,
                      bool m_nofaces = false) {
  // clear();

  if (!m_nofaces) {
    for (typename T2::Finite_faces_iterator it = t2->finite_faces_begin();
         it != t2->finite_faces_end(); ++it) {
      compute_face(it, graphic_buffer, t2, m_color_functor);
    }
  }

  for (typename T2::Finite_edges_iterator it = t2->finite_edges_begin();
       it != t2->finite_edges_end(); ++it) {
    compute_edge<float, T2>(it, graphic_buffer);
  }

  for (typename T2::Finite_vertices_iterator it = t2->finite_vertices_begin();
       it != t2->finite_vertices_end(); ++it) {
    compute_vertex<float, T2>(it, graphic_buffer);
  }
}

} // namespace draw_function_for_t2

template <typename BufferType = float, class T2, class ColorFunctor>
void add_in_graphic_buffer_t2(GraphicBuffer<BufferType> &graphic_buffer,
                              const ColorFunctor &m_color_functor,
                              const T2 *at2 = nullptr, bool m_nofaces = false) {
  if (at2 != nullptr) {
    draw_function_for_t2::compute_elements(graphic_buffer, at2, m_color_functor,
                                           m_nofaces);
  }
}

// Specialization of draw function.
#define CGAL_T2_TYPE CGAL::Triangulation_2<Gt, Tds>

template <class Gt, class Tds,
          class ColorFunctor =
              GenericFunctor<CGAL_T2_TYPE, typename CGAL_T2_TYPE::Vertex_handle,
                             typename CGAL_T2_TYPE::Finite_edges_iterator,
                             typename CGAL_T2_TYPE::Finite_faces_iterator>>
void draw(const CGAL_T2_TYPE &at2,
          const char *title = "Triangulation_2 Basic Viewer",
          const ColorFunctor &color_functor = ColorFunctor(),
          bool nofill = false) {

  GraphicBuffer<float> buffer;
  // add_in_graphic_buffer_t3(buffer, color_functor, &at3, false);
  add_in_graphic_buffer_t2(buffer, DefaultColorFunctorT2(), &at2, false);
  draw_buffer(buffer);
}

#undef CGAL_T2_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_T2_H
