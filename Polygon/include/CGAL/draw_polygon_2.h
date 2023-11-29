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

#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Polygon_2.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
\ingroup PkgDrawPolygon2

opens a new window and draws `ap`, an instance of the `CGAL::Polygon_2` class. A call to this function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam P an instance of the `CGAL::Polygon_2` class.
\param ap the polygon to draw.

*/
template<class P>
void draw(const P& ap);

} /* namespace CGAL */

#endif

namespace CGAL {

namespace draw_function_for_p2 {

template <class P2, class GSOptions>
void compute_elements(const P2& p2,
                      CGAL::Graphics_scene &graphics_scene,
                      const GSOptions& gso)
{
  if (p2.is_empty())
    return;

  typename P2::Point_2 prev=p2.vertex(p2.size()-1);

  if (gso.are_faces_enabled())
  {
    if(gso.colored_face(p2, nullptr))
    { graphics_scene.face_begin(gso.face_color(p2, nullptr)); }
    else
    { graphics_scene.face_begin(); }
  }

  for (typename P2::Vertex_const_iterator i=p2.vertices_begin();
       i!=p2.vertices_end(); ++i)
  {
    if(gso.are_vertices_enabled() &&
       gso.draw_vertex(p2, i))
    { // Add vertex
      if(gso.colored_vertex(p2, i))
      { graphics_scene.add_point(*i, gso.vertex_color(p2, i)); }
      else
      { graphics_scene.add_point(*i); }
    }

    if(gso.are_edges_enabled() &&
       gso.draw_edge(p2, i))
    { // Add edge with previous point
      if(gso.colored_edge(p2, i))
      { graphics_scene.add_segment(prev, *i, gso.edge_color(p2, i)); }
      else
      { graphics_scene.add_segment(prev, *i); }
    }

    if(gso.are_faces_enabled())
    { graphics_scene.add_point_in_face(*i); } // Add point in face

    prev = *i;
  }

  if (gso.are_faces_enabled())
  { graphics_scene.face_end(); }
}

} // namespace draw_function_for_p2

#define CGAL_P2_TYPE CGAL::Polygon_2<T, C>

// Specializations of add_to_graphics_scene function
template<class T, class C, class GSOptions>
void add_to_graphics_scene(const CGAL_P2_TYPE& ap2,
                           CGAL::Graphics_scene& graphics_scene,
                           const GSOptions& gso)
{ draw_function_for_p2::compute_elements(ap2, graphics_scene, gso); }

template<class T, class C>
void add_to_graphics_scene(const CGAL_P2_TYPE& ap2,
                           CGAL::Graphics_scene &graphics_scene)
{
  CGAL::Graphics_scene_options<CGAL_P2_TYPE,
                        typename CGAL_P2_TYPE::Vertex_const_iterator,
                        typename CGAL_P2_TYPE::Vertex_const_iterator,
                        void*> gso;
  draw_function_for_p2::compute_elements(ap2, graphics_scene, gso);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function.
template <class T, class C>
void draw(const CGAL_P2_TYPE &ap2,
          const char *title="Polygon_2 Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(ap2, buffer);
  draw_graphics_scene(buffer, title);
}

template <class T, class C, class GSOptions>
void draw(const CGAL_P2_TYPE &ap2,
          const GSOptions& gso,
          const char *title="Polygon_2 Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(ap2, buffer, gso);
  draw_graphics_scene(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

} // End namespace CGAL

#endif // CGAL_DRAW_POLYGON_2_H
