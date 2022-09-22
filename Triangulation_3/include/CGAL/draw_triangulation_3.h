// Copyright (c) 2018  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_DRAW_T3_H
#define CGAL_DRAW_T3_H

#include <CGAL/Generic_functors.h>
#include <CGAL/GraphicBuffer.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/license/Triangulation_3.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_3.h>

namespace CGAL {

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorT3 {
  template <typename T3>
  static CGAL::IO::Color run(const T3 &,
                             const typename T3::Finite_facets_iterator *fh) {
    if (fh == nullptr)                       // use to get the mono color
      return CGAL::IO::Color(100, 125, 200); // R G B between 0-255

    CGAL::Random random((unsigned int)((std::size_t)(&*((*fh)->first)) +
                                       (std::size_t)((*fh)->second)));
    return get_random_color(random);
  }
};

namespace draw_function_for_t3 {

template <typename BufferType = float, class T3, class ColorFunctor>
void compute_face(typename T3::Finite_facets_iterator fh,
                  GraphicBuffer<BufferType> &graphic_buffer, const T3 *t3,
                  const ColorFunctor &m_fcolor) {

  // TODO: change it to m_fcolor.face_color(*t3, fh)
  CGAL::IO::Color c = m_fcolor.run(*t3, &fh);

  graphic_buffer.face_begin(c);

  graphic_buffer.add_point_in_face(
      fh->first->vertex((fh->second + 1) % 4)->point());
  graphic_buffer.add_point_in_face(
      fh->first->vertex((fh->second + 2) % 4)->point());
  graphic_buffer.add_point_in_face(
      fh->first->vertex((fh->second + 3) % 4)->point());

  graphic_buffer.face_end();
}

template <typename BufferType = float, class T3>
void compute_edge(typename T3::Finite_edges_iterator eh,
                  GraphicBuffer<BufferType> &graphic_buffer) {
  graphic_buffer.add_segment(eh->first->vertex(eh->second)->point(),
                             eh->first->vertex(eh->third)->point());
}

template <typename BufferType = float, class T3>
void compute_vertex(typename T3::Vertex_handle vh,
                    GraphicBuffer<BufferType> &graphic_buffer) {
  graphic_buffer.add_point(vh->point());
}

template <typename BufferType = float, class T3, class ColorFunctor>
void compute_elements(GraphicBuffer<BufferType> &graphic_buffer, const T3 *t3,
                      const ColorFunctor &m_color_functor,
                      bool m_nofaces = false) {

  if (!m_nofaces) {
    for (typename T3::Finite_facets_iterator it = t3->finite_facets_begin();
         it != t3->finite_facets_end(); ++it) {
      compute_face(it, graphic_buffer, t3, m_color_functor);
    }
  }

  for (typename T3::Finite_edges_iterator it = t3->finite_edges_begin();
       it != t3->finite_edges_end(); ++it) {
    compute_edge<float, T3>(it, graphic_buffer);
  }

  for (typename T3::Finite_vertices_iterator it = t3->finite_vertices_begin();
       it != t3->finite_vertices_end(); ++it) {
    compute_vertex<float, T3>(it, graphic_buffer);
  }
}

} // namespace draw_function_for_t3

template <typename BufferType = float, class T3, class ColorFunctor>
void add_in_graphic_buffer_t3(GraphicBuffer<BufferType> &graphic_buffer,
                              const ColorFunctor &m_color_functor,
                              const T3 *at3 = nullptr, bool m_nofaces = false) {
  if (at3 != nullptr) {
    draw_function_for_t3::compute_elements(graphic_buffer, at3, m_color_functor,
                                           m_nofaces);
  }
}

// Specialization of draw function.
#define CGAL_T3_TYPE CGAL::Triangulation_3<Gt, Tds, Lock_data_structure>

template <class Gt, class Tds, class Lock_data_structure,
          class ColorFunctor =
              GenericFunctor<CGAL_T3_TYPE, typename CGAL_T3_TYPE::Vertex_handle,
                             typename CGAL_T3_TYPE::Finite_edges_iterator,
                             typename CGAL_T3_TYPE::Finite_facets_iterator>>
void draw(const CGAL_T3_TYPE &at3, const char *title = "T3 Basic Viewer",
          const ColorFunctor &color_functor = ColorFunctor(),
          bool nofill = false) {

  GraphicBuffer<float> buffer;
  // add_in_graphic_buffer_t3(buffer, color_functor, &at3, false);
  add_in_graphic_buffer_t3(buffer, DefaultColorFunctorT3(), &at3, false);

  draw_buffer(buffer);
}

#undef CGAL_T3_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_T3_H
