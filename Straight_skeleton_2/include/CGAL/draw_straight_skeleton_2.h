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

#include <CGAL/license/Straight_skeleton_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Graphic_storage.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Straight_skeleton_2.h>
#include <sstream>

namespace CGAL {

namespace draw_function_for_ss2 {

template <typename BufferType=float, class SS2, class DrawingFunctor>
void compute_edge(const SS2& ss2, typename SS2::Halfedge_const_handle eh,
                  CGAL::Graphic_storage<BufferType>& graphic_buffer,
                  const DrawingFunctor& drawing_functor)
{
  if (!drawing_functor.draw_edge(ss2, eh))
  { return; }

  if (drawing_functor.colored_edge(ss2, eh))
  {
    graphic_buffer.add_segment(eh->opposite()->vertex()->point(),
                               eh->vertex()->point(),
                               drawing_functor.edge_color(ss2, eh));
  }
  else
  {
    graphic_buffer.add_segment(eh->opposite()->vertex()->point(),
                               eh->vertex()->point());
  }
}

template<typename BufferType=float, class SS2, class DrawingFunctor>
void print_halfedge_labels(const SS2& ss2,
                           typename SS2::Halfedge_const_handle h,
                           CGAL::Graphic_storage<BufferType>& graphic_buffer,
                           const DrawingFunctor& drawing_functor)
{
  // TODO? a functor different from draw_edge to allow to show only some labels ??
  if(!drawing_functor.draw_edge(ss2, h))
  { return; }

  std::stringstream label;
  label << "H" << h->id() << " (V" << h->vertex()->id() << ") ";
  label << "H" << h->opposite()->id() << " (V" << h->opposite()->vertex()->id()
        << ") ";

  graphic_buffer.add_text(
      CGAL::midpoint(h->opposite()->vertex()->point(), h->vertex()->point()),
      label.str());
}

template <typename BufferType = float, class SS2, class DrawingFunctor>
void compute_vertex(const SS2& ss2, typename SS2::Vertex_const_handle vh,
                    CGAL::Graphic_storage<BufferType>& graphic_buffer,
                    const DrawingFunctor& drawing_functor)
{
  if (!drawing_functor.draw_vertex(ss2, vh))
  { return; }

  if (drawing_functor.colored_vertex(ss2, vh))
  {
    graphic_buffer.add_point(vh->point(), drawing_functor.vertex_color(ss2, vh));
  }
  else
  { graphic_buffer.add_point(vh->point()); }
}

template <typename BufferType=float, class SS2, class DrawingFunctor>
void print_vertex_label(const SS2& ss2,
                        typename SS2::Vertex_const_handle vh,
                        CGAL::Graphic_storage<BufferType>& graphic_buffer,
                        const DrawingFunctor& drawing_functor)
{
  // TODO? a functor different from draw_vertex to allow to show only some labels ??
  if (!drawing_functor.draw_vertex(ss2, vh))
  { return; }

  std::stringstream label;
  label << "V" << vh->id() << std::ends;
  graphic_buffer.add_text(vh->point(), label.str());
}

template <typename BufferType=float, class SS2, class DrawingFunctor>
void compute_elements(const SS2& ss2,
                      CGAL::Graphic_storage<BufferType>& graphic_buffer,
                      const DrawingFunctor& drawing_functor)
{
  if (!drawing_functor.are_edges_enabled())
  {
    for (typename SS2::Halfedge_const_iterator it=ss2.halfedges_begin();
         it != ss2.halfedges_end(); ++it)
    {
      if (it->id()<it->opposite()->id())
      {
        compute_edge(ss2, it, graphic_buffer, drawing_functor);
        print_halfedge_labels(ss2, it, graphic_buffer, drawing_functor);
      }
    }
  }

  if (!drawing_functor.are_vertices_enabled())
  {
    for (typename SS2::Vertex_const_iterator it=ss2.vertices_begin();
         it!=ss2.vertices_end(); ++it)
    {
      compute_vertex(ss2, it, graphic_buffer, drawing_functor);
      print_vertex_label(ss2, it, graphic_buffer, drawing_functor);
    }
  }
}

} // namespace draw_function_for_ss2

#define CGAL_SS_TYPE CGAL::Straight_skeleton_2<K>

template <typename BufferType=float, class K, class DrawingFunctor>
void add_in_graphic_buffer(const CGAL_SS_TYPE &ass2,
                           CGAL::Graphic_storage<BufferType>& graphic_buffer,
                           const DrawingFunctor& drawing_functor)
{
  draw_function_for_ss2::compute_elements(ass2, graphic_buffer,
                                          drawing_functor);
}

template <typename BufferType=float, class K>
void add_in_graphic_buffer(const CGAL_SS_TYPE& ass2,
                           CGAL::Graphic_storage<BufferType>& graphic_buffer)
{
  Drawing_functor<CGAL_SS_TYPE,
                  typename CGAL_SS_TYPE::Vertex_const_handle,
                  typename CGAL_SS_TYPE::Halfedge_const_handle,
                  typename CGAL_SS_TYPE::Face_const_handle>
    drawingFunctor;

  drawingFunctor.colored_edge = []
    (const CGAL_SS_TYPE&, typename CGAL_SS_TYPE::Halfedge_const_handle) -> bool
  { return true; };

  drawingFunctor.colored_vertex = []
    (const CGAL_SS_TYPE&, typename CGAL_SS_TYPE::Vertex_const_handle) -> bool
  { return true; };

  drawingFunctor.edge_color = []
    (const CGAL_SS_TYPE&, typename CGAL_SS_TYPE::Halfedge_const_handle eh) -> CGAL::IO::Color
  {
    if (eh->is_bisector())
    { return CGAL::IO::red(); }

    return CGAL::IO::black();
  };

  drawingFunctor.vertex_color = []
    (const CGAL_SS_TYPE&, typename CGAL_SS_TYPE::Vertex_const_handle vh) -> CGAL::IO::Color
  {
    if (vh->is_split())
    { return CGAL::IO::Color(10, 10, 180); } // blue, but not flashy
    else if (vh->has_infinite_time())
    { return CGAL::IO::orange(); }

     return CGAL::IO::Color(10, 180, 10); // green, but not flashy
  };

  add_in_graphic_buffer(ass2, graphic_buffer, drawingFunctor);
}

// Specialization of draw function.
#ifdef CGAL_USE_BASIC_VIEWER

template <class K, class DrawingFunctor>
void draw(const CGAL_SS_TYPE &ass2, const DrawingFunctor &drawingfunctor,
          const char *title="Straight Skeleton Basic Viewer")
{
  CGAL::Graphic_storage<float> buffer;
  add_in_graphic_buffer(ass2, buffer, drawingfunctor);
  draw_buffer(buffer, title);
}

template <class K>
void draw(const CGAL_SS_TYPE &ass2,
          const char *title="Straight Skeleton Basic Viewer")
{
  CGAL::Graphic_storage<float> buffer;
  add_in_graphic_buffer(ass2, buffer);
  draw_buffer(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_SS_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_SS2_H
