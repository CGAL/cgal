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

#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Graphic_storage.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Polygon_2.h>

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

namespace CGAL {

namespace draw_function_for_p2 {

template <typename BufferType=float, class P2, class GSOptions>
void compute_elements(const P2& p2,
                      CGAL::Graphic_storage<BufferType> &graphic_storage,
                      const GSOptions& gs_options)
{
  if (p2.is_empty())
    return;

  typename P2::Point_2 prev=p2.vertex(p2.size()-1);

  if (gs_options.are_faces_enabled())
  {
    if(gs_options.colored_face(p2, nullptr))
    { graphic_storage.face_begin(gs_options.face_color(p2, nullptr)); }
    else
    { graphic_storage.face_begin(); }
  }

  for (typename P2::Vertex_const_iterator i=p2.vertices_begin();
       i!=p2.vertices_end(); ++i)
  {
    if(gs_options.are_vertices_enabled() &&
       gs_options.draw_vertex(p2, i))
    { // Add vertex
      if(gs_options.colored_vertex(p2, i))
      { graphic_storage.add_point(*i, gs_options.vertex_color(p2, i)); }
      else
      { graphic_storage.add_point(*i); }
    }

    if(gs_options.are_edges_enabled() &&
       gs_options.draw_edge(p2, i))
    { // Add edge with previous point
      if(gs_options.colored_vertex(p2, i))
      { graphic_storage.add_segment(prev, *i, gs_options.edge_color(p2, i)); }
      else
      { graphic_storage.add_segment(prev, *i); }
    }

    if(gs_options.are_faces_enabled())
    { graphic_storage.add_point_in_face(*i); } // Add point in face

    prev = *i;
  }

  if (gs_options.are_faces_enabled())
  { graphic_storage.face_end(); }
}

} // namespace draw_function_for_p2

#define CGAL_P2_TYPE CGAL::Polygon_2<T, C>

// Specializations of add_in_graphic_storage function

template<typename BufferType=float, class T, class C, class GSOptions>
void add_in_graphic_storage(const CGAL_P2_TYPE& ap2,
                           CGAL::Graphic_storage<BufferType>& graphic_storage,
                           const GSOptions& gs_options)
{ draw_function_for_p2::compute_elements(ap2, graphic_storage, gs_options); }

template<typename BufferType=float, class T, class C>
void add_in_graphic_storage(const CGAL_P2_TYPE& ap2,
                           CGAL::Graphic_storage<BufferType> &graphic_storage)
{
  CGAL::Graphics_scene_options<CGAL_P2_TYPE,
                        typename CGAL_P2_TYPE::Vertex_const_iterator,
                        typename CGAL_P2_TYPE::Vertex_const_iterator,
                        void*> gs_options;
  draw_function_for_p2::compute_elements(ap2, graphic_storage, gs_options);
}

// Specialization of draw function.

#ifdef CGAL_USE_BASIC_VIEWER

template <class T, class C>
void draw(const CGAL_P2_TYPE &ap2,
          const char *title="Polygon_2 Basic Viewer")
{
  CGAL::Graphic_storage<float> buffer;
  add_in_graphic_storage(ap2, buffer);
  draw_graphic_storage(buffer, title);
}

template <class T, class C, class GSOptions>
void draw(const CGAL_P2_TYPE &ap2,
          const GSOptions& gs_options,
          const char *title="Polygon_2 Basic Viewer")
{
  CGAL::Graphic_storage<float> buffer;
  add_in_graphic_storage(ap2, buffer, gs_options);
  draw_graphic_storage(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

} // End namespace CGAL

#endif // CGAL_DRAW_POLYGON_2_H
