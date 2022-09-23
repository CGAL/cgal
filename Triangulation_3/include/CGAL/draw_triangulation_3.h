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
#include <CGAL/Graphic_buffer.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/license/Triangulation_3.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_3.h>

namespace CGAL {

namespace draw_function_for_t3
{

template <typename BufferType = float, class T3, class DrawingFunctor>
void compute_face(typename T3::Finite_facets_iterator fh,
                  const DrawingFunctor &drawing_functor,
                  CGAL::GraphicBuffer<BufferType> &graphic_buffer, const T3 *t3)
{
  // DONE: TEST draw_face
  if(!drawing_functor.draw_face(*t3, fh))
  {
    return;
  }

  if(drawing_functor.face_color)
  {
    CGAL::IO::Color c = drawing_functor.face_color(*t3, fh);
    graphic_buffer.face_begin(c);
  }
  else
  { graphic_buffer.face_begin(); }

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
                  CGAL::GraphicBuffer<BufferType> &graphic_buffer, const T3 *t3) {
  // DONE: TEST draw_edge + edge color
  if(!drawing_functor.draw_edge(*t3, eh)) {
      return;
  }
  if( && drawing_functor.edge_color)  {

    CGAL::IO::Color c = drawing_functor.edge_color(*t3, eh);
    graphic_buffer.add_segment(eh->first->vertex(eh->second)->point(),
                             eh->first->vertex(eh->third)->point(), c);
  }
  else {
    graphic_buffer.add_segment(eh->first->vertex(eh->second)->point(),
                              eh->first->vertex(eh->third)->point());
  }
}

template <typename BufferType = float, class T3>
void compute_vertex(typename T3::Vertex_handle vh,
                    CGAL::GraphicBuffer<BufferType> &graphic_buffer, const T3 *t3) {

  // DONE TEST draw_vertex  + vertex color
  if(!drawing_functor.draw_vertex(*t3, vh))
  {
    return;
  }

  if(drawing_functor.vertex_color)
  {
    CGAL::IO::Color c = drawing_functor.vertex_color(*t3, vh);
    graphic_buffer.add_point(vh->point(), c);
  }
  else
  { graphic_buffer.add_point(vh->point()); }

}

template <typename BufferType = float, class T3, class DrawingFunctor>
void compute_elements(CGAL::GraphicBuffer<BufferType> &graphic_buffer, const T3 *t3,
                      const DrawingFunctor &drawing_functor,
                      bool m_nofaces = false) {

  if (!m_nofaces) {
    for (typename T3::Finite_facets_iterator it = t3->finite_facets_begin();
         it != t3->finite_facets_end(); ++it)
    {
      compute_face(it, drawing_functor, graphic_buffer, t3);
    }
  }

  for (typename T3::Finite_edges_iterator it = t3->finite_edges_begin();
       it != t3->finite_edges_end(); ++it)
  {
    compute_edge<float, T3>(it, graphic_buffer);
  }

  for (typename T3::Finite_vertices_iterator it = t3->finite_vertices_begin();
       it != t3->finite_vertices_end(); ++it)
  {
    compute_vertex<float, T3>(it, graphic_buffer);
  }
}

} // namespace draw_function_for_t3

template <typename BufferType = float, class T3, class DrawingFunctor>
void add_in_graphic_buffer_t3(CGAL::GraphicBuffer<BufferType> &graphic_buffer,
                              const DrawingFunctor &drawing_functor,
                              const T3 *at3 = nullptr)
{
  if (at3 != nullptr)
  {
    draw_function_for_t3::compute_elements(graphic_buffer, at3, drawing_functor);
  }
}

// Specialization of draw function.
#define CGAL_T3_TYPE CGAL::Triangulation_3<Gt, Tds, Lock_data_structure>

template<class Gt, class Tds, class Lock_data_structure, class DrawingFunctor>
void draw(const CGAL_T3_TYPE &at3, const DrawingFunctor &drawingfunctor,
          const char *title = "T3 Basic Viewer")
{
  CGAL::GraphicBuffer<float> buffer;
  add_in_graphic_buffer_t3(buffer, drawingfunctor, &at3);
  draw_buffer(buffer);
}

template <class Gt, class Tds, class Lock_data_structure>
void draw(const CGAL_T3_TYPE &at3, const char *title = "T3 Basic Viewer")
{
  CGAL::GraphicBuffer<float> buffer;
  CGAL::DrawingFunctor<CGAL_T3_TYPE,
                       typename CGAL_T3_TYPE::Vertex_handle,
                       typename CGAL_T3_TYPE::Finite_edges_iterator,
                       typename CGAL_T3_TYPE::Finite_facets_iterator>
    drawingfunctor;

  drawingfunctor.face_color =
    [](const CGAL_T3_TYPE &at3, const typename CGAL_T3_TYPE::Finite_facets_iterator fh)
    -> CGAL::IO::Color
    {
     if (fh==at3.finite_facets_end())         // use to get the mono color
       return CGAL::IO::Color(100, 125, 200); // R G B between 0-255

     CGAL::Random random((unsigned int)((std::size_t)(&*(fh->first)) +
                                        (std::size_t)(fh->second)));
     
     return get_random_color(random);
     };

  draw(at3, drawingfunctor, title);
}

#undef CGAL_T3_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_T3_H
