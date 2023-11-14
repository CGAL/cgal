// Copyright (c) 2022 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s):   Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//              Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_GRAPHICS_SCENE_OPTIONS_H
#define CGAL_GRAPHICS_SCENE_OPTIONS_H

#include <CGAL/license/GraphicsView.h>
#include <CGAL/IO/Color.h>
#include <functional>

namespace CGAL {

template <typename DS,
          typename VertexHandle,
          typename EdgeHandle,
          typename FaceHandle,
          typename VolumeHandle=void>
struct Graphics_scene_options;

// Drawing functor for a 2D combinatorial data structure
// (with vertices, edges and faces)
template <typename DS,
          typename VertexHandle,
          typename EdgeHandle,
          typename FaceHandle>
struct Graphics_scene_options<DS, VertexHandle, EdgeHandle, FaceHandle, void>
{

  typedef VertexHandle vertex_descriptor;
  typedef EdgeHandle edge_descriptor;
  typedef FaceHandle face_descriptor;

  Graphics_scene_options(): m_enabled_vertices(true),
                            m_enabled_edges(true),
                            m_enabled_faces(true)
  {
    draw_vertex=[](const DS &, VertexHandle)->bool { return true; };
    draw_edge=[](const DS &, EdgeHandle)->bool { return true; };
    draw_face=[](const DS &, FaceHandle)->bool { return true; };

    colored_vertex=[](const DS &, VertexHandle)->bool { return false; };
    colored_edge=[](const DS &, EdgeHandle)->bool { return false; };
    colored_face=[](const DS &, FaceHandle)->bool { return false; };

    face_wireframe=[](const DS &, FaceHandle)->bool { return false; };
  }

  // The seven following functions should not be null
  std::function<bool(const DS &, VertexHandle)> draw_vertex;
  std::function<bool(const DS &, EdgeHandle)>   draw_edge;
  std::function<bool(const DS &, FaceHandle)>   draw_face;

  std::function<bool(const DS &, VertexHandle)> colored_vertex;
  std::function<bool(const DS &, EdgeHandle)>   colored_edge;
  std::function<bool(const DS &, FaceHandle)>   colored_face;

  std::function<bool(const DS &, FaceHandle)> face_wireframe;

  // These functions must be non null if the corresponding colored_XXX function
  // returns true.
  std::function<CGAL::IO::Color(const DS &, VertexHandle)> vertex_color;
  std::function<CGAL::IO::Color(const DS &, EdgeHandle)>   edge_color;
  std::function<CGAL::IO::Color(const DS &, FaceHandle)>   face_color;

  void disable_vertices() { m_enabled_vertices=false; }
  void enable_vertices() { m_enabled_vertices=true; }
  bool are_vertices_enabled() const { return m_enabled_vertices; }

  void disable_edges() { m_enabled_edges=false; }
  void enable_edges() { m_enabled_edges=true; }
  bool are_edges_enabled() const { return m_enabled_edges; }

  void disable_faces() { m_enabled_faces=false; }
  void enable_faces() { m_enabled_faces=true; }
  bool are_faces_enabled() const { return m_enabled_faces; }

protected:
  bool m_enabled_vertices, m_enabled_edges, m_enabled_faces;
};

// Drawing functor for a 3D combinatorial data structure
// (with vertices, edges, faces and volumes)
template <typename DS,
          typename VertexHandle,
          typename EdgeHandle,
          typename FaceHandle,
          typename VolumeHandle>
struct Graphics_scene_options:
    public Graphics_scene_options<DS, VertexHandle, EdgeHandle, FaceHandle>
{
  typedef VertexHandle vertex_descriptor;
  typedef EdgeHandle edge_descriptor;
  typedef FaceHandle face_descriptor;
  typedef VolumeHandle volume_descriptor;

  Graphics_scene_options() : m_enabled_volumes(true)
  {
    draw_volume=[](const DS &, VolumeHandle)->bool { return true; };
    colored_volume=[](const DS &, VolumeHandle)->bool { return false; };
    volume_wireframe=[](const DS &, VolumeHandle)->bool { return false; };
  }

  std::function<bool(const DS &, VolumeHandle)>            draw_volume;
  std::function<bool(const DS &, VolumeHandle)>            colored_volume;
  std::function<bool(const DS &, VolumeHandle)>            volume_wireframe;
  std::function<CGAL::IO::Color(const DS &, VolumeHandle)> volume_color;

  void disable_volumes() { m_enabled_volumes=false; }
  void enable_volumes() { m_enabled_volumes=true; }
  bool are_volumes_enabled() const { return m_enabled_volumes; }

protected:
  bool m_enabled_volumes;
};

} // End namespace CGAL

#endif // CGAL_GRAPHICS_SCENE_OPTIONS_H
