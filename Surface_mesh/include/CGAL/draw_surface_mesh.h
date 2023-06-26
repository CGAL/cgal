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

Open a new window and draw `asm`, an instance of the `CGAL::Surface_mesh` class. The function is blocking, that is the program continues as soon as the user closes the window. This function requires `CGAL_Qt5`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt5` and add the definition `CGAL_USE_BASIC_VIEWER`.
\tparam SM an instance of the `CGAL::Surface_mesh` class.
\param asm the surface mesh to draw.

*/
template<class SM>
void draw(const SM& asm);

#else // DOXYGEN_RUNNING

#include <CGAL/license/Surface_mesh.h>
#include <CGAL/Graphic_storage.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_face_graph.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

namespace CGAL {

// Check if there are any color maps that could be used, random otherwise
template <typename K>
struct Surface_mesh_basic_viewer_color_map
  : public Drawing_functor<Surface_mesh<K>,
                           typename boost::graph_traits<::CGAL::Surface_mesh<K>>::vertex_descriptor,
                           typename boost::graph_traits<::CGAL::Surface_mesh<K>>::edge_descriptor,
                           typename boost::graph_traits<::CGAL::Surface_mesh<K>>::face_descriptor>
{
  using SM = ::CGAL::Surface_mesh<K>;
  using vertex_descriptor = typename boost::graph_traits<SM>::vertex_descriptor;
  using edge_descriptor = typename boost::graph_traits<SM>::edge_descriptor;
  using face_descriptor = typename boost::graph_traits<SM>::face_descriptor;

  using Base=Drawing_functor<SM, vertex_descriptor, edge_descriptor, face_descriptor>;

  using vertex_colors = typename SM::template Property_map<vertex_descriptor, CGAL::IO::Color>;
  using edge_colors = typename SM::template Property_map<edge_descriptor, CGAL::IO::Color>;
  using face_colors = typename SM::template Property_map<face_descriptor, CGAL::IO::Color>;

  Surface_mesh_basic_viewer_color_map(const SM& amesh)
  {
    bool found = false;
    std::tie(vcolors, found) = amesh.template property_map<vertex_descriptor, CGAL::IO::Color>("v:color");
    std::tie(ecolors, found) = amesh.template property_map<edge_descriptor, CGAL::IO::Color>("e:color");
    std::tie(fcolors, found) = amesh.template property_map<face_descriptor, CGAL::IO::Color>("f:color");
    CGAL_USE(found);
  }

  CGAL::IO::Color operator()(const Surface_mesh<K>& amesh,
                             const vertex_descriptor v) const
  {
    return vcolors ? get(vcolors, v) : Base::operator()(amesh, v);
  }

  CGAL::IO::Color operator()(const Surface_mesh<K>& amesh,
                             const edge_descriptor e) const
  {
    return ecolors ? get(ecolors, e) : Base::operator()(amesh, e);
  }

  CGAL::IO::Color operator()(const Surface_mesh<K>& amesh,
                             const face_descriptor f) const
  {
    return fcolors ? get(fcolors, f) : Base::operator()(amesh, f);
  }

private:
  vertex_colors vcolors;
  edge_colors ecolors;
  face_colors fcolors;
};

template<class K, typename BufferType=float,  class DrawingFunctor>
void add_in_graphic_storage(const Surface_mesh<K>& amesh,
                           CGAL::Graphic_storage<BufferType> &graphic_storage,
                           const DrawingFunctor &drawing_functor)
{ add_in_graphic_storage_for_fg(amesh, graphic_storage, drawing_functor); }

template<class K, typename BufferType=float>
void add_in_graphic_storage(const Surface_mesh<K>& amesh,
                           CGAL::Graphic_storage<BufferType> &graphic_storage)
{ add_in_graphic_storage_for_fg(amesh, graphic_storage,
                                Surface_mesh_basic_viewer_color_map<K>(amesh)); }

#ifdef CGAL_USE_BASIC_VIEWER

  // Specialization of draw function.
template<class K, typename BufferType=float>
void draw(const Surface_mesh<K>& amesh,
          const char* title="Surface_mesh Basic Viewer")
{
  CGAL::Graphic_storage<BufferType> buffer;
  add_in_graphic_storage_for_fg(amesh, buffer);
  draw_graphic_storage(buffer, title);
}

template<class K, typename BufferType=float, class DrawingFunctor>
void draw(const Surface_mesh<K>& amesh,
          const DrawingFunctor &drawing_functor,
          const char* title="Surface_mesh Basic Viewer")
{
  CGAL::Graphic_storage<BufferType> buffer;
  add_in_graphic_storage_for_fg(amesh, buffer, drawing_functor);
  draw_graphic_storage(buffer, title);
}

#endif // CGAL_USE_BASIC_VIEWER

} // End namespace CGAL

#endif // DOXYGEN_RUNNING

#endif // CGAL_DRAW_SURFACE_MESH_H
