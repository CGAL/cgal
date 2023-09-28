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

#include <CGAL/license/Triangulation_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Graphic_storage.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Random.h>

namespace CGAL {

namespace draw_function_for_t2 {

template <typename BufferType=float, class T2, class GSOptions>
void compute_face(const T2& t2,
                  typename T2::Finite_faces_iterator fh,
                  CGAL::Graphic_storage<BufferType>& graphic_storage,
                  const GSOptions& gs_options)
{
  if (!gs_options.draw_face(t2, fh))
  { return; }

  if (gs_options.colored_face(t2, fh))
  { graphic_storage.face_begin(gs_options.face_color(t2, fh)); }
  else
  { graphic_storage.face_begin(); }

  graphic_storage.add_point_in_face(fh->vertex(0)->point());
  graphic_storage.add_point_in_face(fh->vertex(1)->point());
  graphic_storage.add_point_in_face(fh->vertex(2)->point());

  graphic_storage.face_end();
}

template <typename BufferType=float, class T2, class GSOptions>
void compute_edge(const T2& t2, typename T2::Finite_edges_iterator eh,
                  CGAL::Graphic_storage<BufferType>& graphic_storage,
                  const GSOptions& gs_options)
{
  if (!gs_options.draw_edge(t2, eh))
  { return; }

  if (gs_options.colored_edge(t2, eh))
  {
    graphic_storage.add_segment
      (eh->first->vertex(eh->first->cw(eh->second))->point(),
       eh->first->vertex(eh->first->ccw(eh->second))->point(),
       gs_options.edge_color(t2, eh));
  }
  else
  {
    graphic_storage.add_segment
      (eh->first->vertex(eh->first->cw(eh->second))->point(),
       eh->first->vertex(eh->first->ccw(eh->second))->point());
  }
}

template <typename BufferType=float, class T2, class GSOptions>
void compute_vertex(const T2& t2, typename T2::Vertex_handle vh,
                    CGAL::Graphic_storage<BufferType>& graphic_storage,
                    const GSOptions& gs_options)
{
  if (!gs_options.draw_vertex(t2, vh))
  { return; }

  if (gs_options.colored_vertex(t2, vh))
  {
    graphic_storage.add_point(vh->point(), gs_options.vertex_color(t2, vh));
  }
  else
  {
    graphic_storage.add_point(vh->point());
  }
}

template <typename BufferType=float, class T2, class GSOptions>
void compute_elements(const T2& t2,
                      CGAL::Graphic_storage<BufferType>& graphic_storage,
                      const GSOptions& gs_options)
{
  if (gs_options.are_faces_enabled())
  {
    for (typename T2::Finite_faces_iterator it=t2.finite_faces_begin();
         it!=t2.finite_faces_end(); ++it)
    { compute_face(t2, it, graphic_storage, gs_options); }
  }

  if (gs_options.are_edges_enabled())
  {
    for (typename T2::Finite_edges_iterator it=t2.finite_edges_begin();
         it!=t2.finite_edges_end(); ++it)
    { compute_edge(t2, it, graphic_storage, gs_options); }
  }

  if (gs_options.are_vertices_enabled())
  {
    for (typename T2::Finite_vertices_iterator it=t2.finite_vertices_begin();
         it!=t2.finite_vertices_end(); ++it)
    { compute_vertex(t2, it, graphic_storage, gs_options); }
  }
}

} // namespace draw_function_for_t2

#define CGAL_T2_TYPE CGAL::Triangulation_2<Gt, Tds>

template <class Gt, class Tds, typename BufferType=float, class GSOptions>
void add_in_graphic_storage(const CGAL_T2_TYPE& at2,
                            CGAL::Graphic_storage<BufferType>& graphic_storage,
                            const GSOptions& gs_options)
{
  draw_function_for_t2::compute_elements(at2, graphic_storage, gs_options);
}

template <class Gt, class Tds, typename BufferType=float>
void add_in_graphic_storage(const CGAL_T2_TYPE& at2,
                            CGAL::Graphic_storage<BufferType>& graphic_storage)
{
  Graphics_scene_options<CGAL_T2_TYPE,
                  typename CGAL_T2_TYPE::Vertex_handle,
                  typename CGAL_T2_TYPE::Finite_edges_iterator,
                  typename CGAL_T2_TYPE::Finite_faces_iterator>
      drawingFunctor;

  drawingFunctor.colored_face =
      [](const CGAL_T2_TYPE&, const typename CGAL_T2_TYPE::Finite_faces_iterator) -> bool
      { return true; };

  drawingFunctor.face_color =
      [](const CGAL_T2_TYPE&, const typename CGAL_T2_TYPE::Finite_faces_iterator fh) -> CGAL::IO::Color
      {
        CGAL::Random random((unsigned int)(std::size_t)(&*fh));
        return get_random_color(random);
      };

  add_in_graphic_storage(at2, graphic_storage, drawingFunctor);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function.
template <class Gt, class Tds, class GSOptions>
void draw(const CGAL_T2_TYPE &at2, const GSOptions &gs_options,
          const char *title="Triangulation_2 Basic Viewer")
{
  CGAL::Graphic_storage<float> buffer;
  add_in_graphic_storage(at2, buffer, gs_options);
  draw_graphic_storage(buffer, title);
}

template <class Gt, class Tds>
void draw(const CGAL_T2_TYPE& at2,
          const char *title="Triangulation_2 Basic Viewer")
{
  CGAL::Graphic_storage<float> buffer;
  add_in_graphic_storage(at2, buffer);
  draw_graphic_storage(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_T2_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_T2_H
