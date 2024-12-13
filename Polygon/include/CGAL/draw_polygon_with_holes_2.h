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

#ifndef CGAL_DRAW_POLYGON_WITH_HOLES_2_H
#define CGAL_DRAW_POLYGON_WITH_HOLES_2_H

#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Polygon_with_holes_2.h>

#ifdef DOXYGEN_RUNNING
namespace CGAL {

/*!
  \ingroup PkgDrawPolygonWithHoles2

opens a new window and draws a 2D polygon with holes. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam PH which must be an instantiation of a `CGAL::Polygon_with_holes_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param ph the polygon with holes to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class T, class C, class GSOptions>

 void CGAL::draw(const CGAL::Polygon_with_holes_2<T, C>& ph, const GSOptions& gso);
</code>
\cgalAdvancedEnd
 */
template <typename PH, class GSOptions>
void draw(const PH& ph, const GSOptions& gso);

/*!
\ingroup PkgDrawPolygonWithHoles2

A shortcut to `CGAL::draw(ph, Graphics_scene_options{})`.
*/
  template<class PH>
  void draw(const PH& ph);

/*!
\ingroup PkgDrawPolygonWithHoles2

adds the vertices, edges and faces of `ph` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam PH which must be an instantiation of a `CGAL::Polygon_with_holes_2<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param ph the polygon with holes to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class T, class C, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Polygon_with_holes_2<T, C>& ph, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class PH, class GSOptions>
void add_to_graphics_scene(const PH& ph,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawPolygonWithHoles2

A shortcut to `CGAL::add_to_graphics_scene(ph, gs, Graphics_scene_options{})`.
*/
template<class PH>
void add_to_graphics_scene(const PH& ph,
                           CGAL::Graphics_scene& gs);

} /* namespace CGAL */

#endif

namespace CGAL
{

namespace draw_function_for_ph2_with_holes {

template <class P2, class GSOptions>
void compute_one_loop_elements(const P2& ap2,
                               const typename P2::General_polygon_2& aloop,
                               Graphics_scene &graphics_scene,
                               bool hole,
                               const GSOptions& gs_options)
{
  if (hole && gs_options.are_faces_enabled())
  { graphics_scene.add_point_in_face(aloop.vertex(aloop.size()-1)); }

  typename P2::General_polygon_2::Vertex_const_iterator prev;
  for(typename P2::General_polygon_2::Vertex_const_iterator i=aloop.vertices_begin();
      i!=aloop.vertices_end(); ++i)
  {
    if(gs_options.are_vertices_enabled() &&
       gs_options.draw_vertex(ap2, i))
    { // Add vertex
      if(gs_options.colored_vertex(ap2, i))
      { graphics_scene.add_point(*i, gs_options.vertex_color(ap2, i)); }
      else
      { graphics_scene.add_point(*i); }
    }

    if(i!=aloop.vertices_begin() &&
       gs_options.are_edges_enabled() && gs_options.draw_edge(ap2, prev))
    { // Add segment with previous point
      if(gs_options.colored_edge(ap2, prev))
      { graphics_scene.add_segment(*prev, *i, gs_options.edge_color(ap2, prev)); }
      else
      { graphics_scene.add_segment(*prev, *i); }
    }

    if(gs_options.are_faces_enabled())
    { graphics_scene.add_point_in_face(*i); } // Add point in face

    prev=i;
  }

  // Add the last segment between the last point and the first one
  if(gs_options.are_edges_enabled() &&
     gs_options.draw_edge(ap2, aloop.vertices_begin()))
  {
    if(gs_options.colored_edge(ap2, prev))
    { graphics_scene.add_segment(*prev, *(aloop.vertices_begin()),
                                 gs_options.edge_color(ap2, prev)); }
    else
    { graphics_scene.add_segment(*prev, *(aloop.vertices_begin())); }
  }
}

template <class P2, class GSOptions>
void compute_elements(const P2& p2, Graphics_scene &graphics_scene,
                      const GSOptions& gs_options)
{
  if (p2.outer_boundary().is_empty()) return;

  if (gs_options.are_faces_enabled())
  {
    if(gs_options.colored_face(p2, nullptr))
    { graphics_scene.face_begin(gs_options.face_color(p2, nullptr)); }
    else
    { graphics_scene.face_begin(); }
  }

  compute_one_loop_elements<P2>(p2, p2.outer_boundary(), graphics_scene,
                                false, gs_options);

  for (typename P2::Hole_const_iterator it=p2.holes_begin(); it!=p2.holes_end(); ++it)
  {
    compute_one_loop_elements<P2>(p2, *it, graphics_scene,
                                  true, gs_options);
    if (gs_options.are_faces_enabled())
    { graphics_scene.add_point_in_face(p2.outer_boundary().vertex
                                       (p2.outer_boundary().size()-1));
    }
  }

 if (gs_options.are_faces_enabled())
 { graphics_scene.face_end(); }
}

} // draw_function_for_ph2

#define CGAL_P2_WITH_HOLES_TYPE CGAL::Polygon_with_holes_2<T, C>

template <class T, class C, class GSOptions>
void add_to_graphics_scene(const CGAL_P2_WITH_HOLES_TYPE& p2,
                           CGAL::Graphics_scene& graphics_scene,
                           const GSOptions &gs_options)
{
  draw_function_for_ph2_with_holes::compute_elements(p2, graphics_scene,
                                                     gs_options);
}

template <class T, class C>
void add_to_graphics_scene(const CGAL_P2_WITH_HOLES_TYPE& p2,
                           CGAL::Graphics_scene& graphics_scene)
{
  Graphics_scene_options<CGAL_P2_WITH_HOLES_TYPE,
                  typename CGAL_P2_WITH_HOLES_TYPE::General_polygon_2::Vertex_const_iterator,
                  typename CGAL_P2_WITH_HOLES_TYPE::General_polygon_2::Vertex_const_iterator,
                  void*> gs_options;

  add_to_graphics_scene(p2, graphics_scene, gs_options);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function.
template<class T, class C, class GSOptions>
void draw(const CGAL_P2_WITH_HOLES_TYPE& ap2, const GSOptions &gs_options,
          const char* title="Polygon with Holes Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(ap2, buffer, gs_options);
  draw_graphics_scene(buffer, title);
}

template<class T, class C>
void draw(const CGAL_P2_WITH_HOLES_TYPE& ap2,
          const char* title="Polygon with Holes Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(ap2, buffer);
  draw_graphics_scene(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_P2_WITH_HOLES_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_POLYGON_WITH_HOLES_2_H
