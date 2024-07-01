// Copyright (c) 2018-2020 GeometryFactory (France)
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

#ifndef CGAL_DRAW_SURFACE_MESH_H
#define CGAL_DRAW_SURFACE_MESH_H

#ifdef DOXYGEN_RUNNING

/*!
\ingroup PkgDrawSurfaceMesh

opens a new window and draws a surface mesh. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam SM which must be an instantiation of a `CGAL::Surface_mesh<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param sm the surface mesh to draw.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class K, class GSOptions>

 void CGAL::draw(const CGAL::Surface_mesh<K>& sm, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class SM, class GSOptions>
void draw(const SM& sm, const GSOptions& gso);

/*!
\ingroup PkgDrawSurfaceMesh

A shortcut to `CGAL::draw(sm, Graphics_scene_options{})`.
*/
  template<class SM>
  void draw(const SM& sm);

/*!
\ingroup PkgDrawSurfaceMesh

adds the vertices, edges and faces of `sm` into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam SM which must be an instantiation of a `CGAL::Surface_mesh<...>`.
\tparam GSOptions a model of `GraphicsSceneOptions` concept.

\param sm the surface mesh to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.

\cgalAdvancedBegin
The real declaration of this function template is:

<code>
 template<class K, class GSOptions>

 void CGAL::add_to_graphics_scene(const CGAL::Surface_mesh<K>& sm, CGAL::Graphics_scene& gs, const GSOptions& gso);
</code>
\cgalAdvancedEnd
*/
template<class SM, class GSOptions>
void add_to_graphics_scene(const SM& sm,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso);

/*!
\ingroup PkgDrawSurfaceMesh

A shortcut to `CGAL::add_to_graphics_scene(sm, gs, Graphics_scene_options{})`.
*/
template<class SM>
void add_to_graphics_scene(const SM& sm,
                           CGAL::Graphics_scene& gs);

#else // DOXYGEN_RUNNING

#include <CGAL/license/Surface_mesh.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_face_graph.h>
#include <CGAL/Qt/Basic_viewer.h>

namespace CGAL {

// Check if there are any color maps that could be used
template <typename K>
struct Graphics_scene_options_surface_mesh
  : public Graphics_scene_options<Surface_mesh<K>,
                           typename boost::graph_traits<::CGAL::Surface_mesh<K>>::vertex_descriptor,
                           typename boost::graph_traits<::CGAL::Surface_mesh<K>>::edge_descriptor,
                           typename boost::graph_traits<::CGAL::Surface_mesh<K>>::face_descriptor>
{
  using SM = ::CGAL::Surface_mesh<K>;
  using vertex_descriptor = typename boost::graph_traits<SM>::vertex_descriptor;
  using edge_descriptor = typename boost::graph_traits<SM>::edge_descriptor;
  using face_descriptor = typename boost::graph_traits<SM>::face_descriptor;

  Graphics_scene_options_surface_mesh(const SM& amesh)
  {
    auto _vcolors = amesh.template property_map<vertex_descriptor, CGAL::IO::Color>("v:color");
    if(_vcolors.has_value())
    {
      vcolors = _vcolors.value();
      this->colored_vertex=[](const SM &, vertex_descriptor)->bool { return true; };
      this->vertex_color=[this](const SM &, vertex_descriptor v)->CGAL::IO::Color
      { return get(vcolors, v); };
    }
    else
    { this->colored_vertex=[](const SM &, vertex_descriptor)->bool { return false; }; }

    auto _ecolors = amesh.template property_map<edge_descriptor, CGAL::IO::Color>("e:color");
    if(_ecolors.has_value())
    {
      ecolors = _ecolors.value();
      this->colored_edge=[](const SM &, edge_descriptor)->bool { return true; };
      this->edge_color=[this](const SM &, edge_descriptor e)->CGAL::IO::Color
      { return get(ecolors, e); };
    }
    else
    { this->colored_edge=[](const SM &, edge_descriptor)->bool { return false; }; }

    auto _fcolors = amesh.template property_map<face_descriptor, CGAL::IO::Color>("f:color");
    if(_fcolors.has_value())
    {
      fcolors = _fcolors.value();
      this->colored_face=[](const SM &, face_descriptor)->bool { return true; };
      this->face_color=[this](const SM &, face_descriptor f)->CGAL::IO::Color
      { return get(fcolors, f); };
    }
    else
    { this->colored_face=[](const SM &, face_descriptor)->bool { return false; }; }
  }

private:
  typename SM::template Property_map<vertex_descriptor, CGAL::IO::Color> vcolors;
  typename SM::template Property_map<edge_descriptor, CGAL::IO::Color> ecolors;
  typename SM::template Property_map<face_descriptor, CGAL::IO::Color> fcolors;
};

template<class K,  class GSOptions>
void add_to_graphics_scene(const Surface_mesh<K>& amesh,
                           CGAL::Graphics_scene &graphics_scene,
                           const GSOptions &gs_options)
{ add_to_graphics_scene_for_fg(amesh, graphics_scene, gs_options); }

template<class K>
void add_to_graphics_scene(const Surface_mesh<K>& amesh,
                           CGAL::Graphics_scene &graphics_scene)
{ add_to_graphics_scene_for_fg(amesh, graphics_scene,
                               Graphics_scene_options_surface_mesh<K>(amesh)); }

#ifdef CGAL_USE_BASIC_VIEWER

  // Specialization of draw function.
template<class K>
void draw(const Surface_mesh<K>& amesh,
          const char* title="Surface_mesh Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(amesh, buffer);
  draw_graphics_scene(buffer, title);
}

template<class K, class GSOptions>
void draw(const Surface_mesh<K>& amesh,
          const GSOptions &gs_options,
          const char* title="Surface_mesh Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(amesh, buffer, gs_options);
  draw_graphics_scene(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

} // End namespace CGAL

#endif // DOXYGEN_RUNNING

#endif // CGAL_DRAW_SURFACE_MESH_H
