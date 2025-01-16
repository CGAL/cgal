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

opens a new window and draws a 2D polygon. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam P which must be an instantiation of a `CGAL::Polygon_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param p the polygon to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class T, class C, class GSOptions>

 void CGAL::draw(const CGAL::Polygon_2<T, C>& p, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class P, class GSOptions>
void draw(const P& p, const GSOptions& gso);

/*!
\ingroup PkgDrawPolygon2

A shortcut to `CGAL::draw(p, Graphics_scene_options{})`.
*/
  template<class P>
  void draw(const P& p);

/*!
\ingroup PkgDrawPolygon2

adds the vertices, edges and faces of `p` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam P which must be an instantiation of a `CGAL::Polygon_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param p the 2D polygon to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class T, class C, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Polygon_2<T, C>& p, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class P, class GSOptions>
void add_to_graphics_scene(const P& p,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawPolygon2

A shortcut to `CGAL::add_to_graphics_scene(p, gs, Graphics_scene_options{})`.
*/
template<class P>
void add_to_graphics_scene(const P& p,
                           CGAL::Graphics_scene& gs);
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

  if (gso.are_faces_enabled())
  {
    if(gso.colored_face(p2, nullptr))
    { graphics_scene.face_begin(gso.face_color(p2, nullptr)); }
    else
    { graphics_scene.face_begin(); }
  }

  typename P2::Vertex_const_iterator prev=p2.vertices_end(); --prev;
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
       gso.draw_edge(p2, prev))
    { // Add edge with previous point
      if(gso.colored_edge(p2, prev))
      { graphics_scene.add_segment(*prev, *i, gso.edge_color(p2, prev)); }
      else
      { graphics_scene.add_segment(*prev, *i); }
    }

    if(gso.are_faces_enabled())
    { graphics_scene.add_point_in_face(*i); } // Add point in face

    prev=i;
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

#undef CGAL_P2_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_POLYGON_2_H
