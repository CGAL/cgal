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

#ifndef CGAL_DRAW_SS2_H
#define CGAL_DRAW_SS2_H

#include <CGAL/Drawing_functor.h>
#include <CGAL/Graphic_buffer.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/license/Straight_skeleton_2.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Random.h>
#include <CGAL/Straight_skeleton_2.h>

#include <sstream>

namespace CGAL {

namespace draw_function_for_ss2 {

template <typename BufferType = float, class SS2, class DrawingFunctor>
void compute_edge(typename SS2::Halfedge_const_handle eh,
                  const DrawingFunctor &drawing_functor, const SS2 *ss2,
                  CGAL::GraphicBuffer<BufferType> &graphic_buffer) {

  if (!drawing_functor.draw_edge(*ss2, eh)) {
    return;
  }

  if (drawing_functor.colored_edge(*ss2, eh)) {
    if (eh->is_bisector())
      graphic_buffer.add_segment(eh->opposite()->vertex()->point(),
                                 eh->vertex()->point(), CGAL::IO::red());
    else
      graphic_buffer.add_segment(eh->opposite()->vertex()->point(),
                                 eh->vertex()->point(), CGAL::IO::black());
  }
}

template <typename BufferType = float, class SS2>
void print_halfedge_labels(typename SS2::Halfedge_const_handle h,
                           CGAL::GraphicBuffer<BufferType> &graphic_buffer) {
  std::stringstream label;

  label << "H" << h->id() << " (V" << h->vertex()->id() << ") ";
  label << "H" << h->opposite()->id() << " (V" << h->opposite()->vertex()->id()
        << ") ";

  graphic_buffer.add_text(
      CGAL::midpoint(h->opposite()->vertex()->point(), h->vertex()->point()),
      label.str());
}

template <typename BufferType = float, class SS2, class DrawingFunctor>
void compute_vertex(typename SS2::Vertex_const_handle vh,
                    const DrawingFunctor &drawing_functor, const SS2 *ss2,
                    CGAL::GraphicBuffer<BufferType> &graphic_buffer) {

  if (!drawing_functor.draw_vertex(*ss2, vh)) {
    return;
  }
  if (drawing_functor.colored_vertex(*ss2, vh)) {
    if (vh->is_split())
      graphic_buffer.add_point(
          vh->point(), CGAL::IO::Color(10, 10, 180)); // blue, but not flashy
    else if (vh->has_infinite_time())
      graphic_buffer.add_point(vh->point(), CGAL::IO::orange());
    else
      graphic_buffer.add_point(
          vh->point(), CGAL::IO::Color(10, 180, 10)); // green, but not flashy
  }
}

template <typename BufferType = float, class SS2>
void print_vertex_label(typename SS2::Vertex_const_handle vh,
                        CGAL::GraphicBuffer<BufferType> &graphic_buffer) {
  std::stringstream label;
  label << "V" << vh->id() << std::ends;
  graphic_buffer.add_text(vh->point(), label.str());
}

template <typename BufferType = float, class SS2, class DrawingFunctor>
void compute_elements(const SS2 *ss2,
                      CGAL::GraphicBuffer<BufferType> &graphic_buffer,
                      const DrawingFunctor &drawing_functor) {

  for (typename SS2::Halfedge_const_iterator it = ss2->halfedges_begin();
       it != ss2->halfedges_end(); ++it) {
    if (it->id() < it->opposite()->id()) {
      compute_edge(it, drawing_functor, ss2, graphic_buffer);
      print_halfedge_labels<BufferType, SS2>(it, graphic_buffer);
    }
  }
  for (typename SS2::Vertex_const_iterator it = ss2->vertices_begin();
       it != ss2->vertices_end(); ++it) {
    compute_vertex(it, drawing_functor, ss2, graphic_buffer);
    print_vertex_label<BufferType, SS2>(it, graphic_buffer);
  }
}

} // namespace draw_function_for_ss2

template <typename BufferType = float, class SS2, class DrawingFunctor>
void add_in_graphic_buffer_ss2(const SS2 &ass2,
                               CGAL::GraphicBuffer<BufferType> &graphic_buffer,
                               const DrawingFunctor &drawing_functor) {
  draw_function_for_ss2::compute_elements(&ass2, graphic_buffer,
                                          drawing_functor);
}

template <typename BufferType = float, class SS2>
void add_in_graphic_buffer_ss2(
    const SS2 &ass2, CGAL::GraphicBuffer<BufferType> &graphic_buffer) {
  Drawing_functor<SS2, typename SS2::Vertex_const_handle,
                  typename SS2::Halfedge_const_handle,
                  typename SS2::Face_const_handle>
      drawingFunctor;

  add_in_graphic_buffer_ss2(ass2, graphic_buffer, drawingFunctor);
}

// Specialization of draw function.
#define CGAL_SS_TYPE CGAL::Straight_skeleton_2<K>

template <class K, class DrawingFunctor>
void draw(const CGAL_SS_TYPE &ass2, const DrawingFunctor &drawingfunctor,
          const char *title = "Straight Skeleton Basic Viewer") {
  CGAL::GraphicBuffer<float> buffer;
  add_in_graphic_buffer_ss2(ass2, buffer, drawingfunctor);
  draw_buffer(buffer);
}

template <class K>
void draw(const CGAL_SS_TYPE &ass2,
          const char *title = "Straight Skeleton Basic Viewer") {
  CGAL::GraphicBuffer<float> buffer;
  add_in_graphic_buffer_ss2(ass2, buffer);
  draw_buffer(buffer);
}

#undef CGAL_SS_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_SS2_H
