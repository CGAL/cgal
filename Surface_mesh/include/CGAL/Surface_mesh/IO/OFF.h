// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_IO_OFF_H
#define CGAL_SURFACE_MESH_IO_OFF_H

#include <CGAL/license/Surface_mesh.h>

#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>

#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/IO/Color.h>

#include <iostream>
#include <tuple>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

/// \relates Surface_mesh
///
/// Extracts the surface mesh from an input stream in Ascii OFF, COFF, NOFF, CNOFF
/// format and appends it to the surface mesh `sm`.
///
/// The operator reads the point property as well as "v:normal", "v:color", `"v:texcoord"`, and "f:color".
/// If an alternative vertex_point map is given through `np`,
/// then it will be used instead of the default one.
///
/// \pre The data in the stream must represent a two-manifold. If this is not the case
///      the `failbit` of `is` is set and the mesh cleared.
///
template <typename Point, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(std::istream& is,
              Surface_mesh<Point>& sm,
              const CGAL_BGL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                            Mesh;
  typedef typename Mesh::Vertex_index                                    Vertex_index;
  typedef typename Mesh::Face_index                                      Face_index;

  typedef typename Kernel_traits<Point>::Kernel                          K;
  typedef typename K::Vector_3                                           Normal;
  typedef typename K::Point_2                                            Texture;
  typedef CGAL::Color                                                    Color;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typename CGAL::GetVertexPointMap<Mesh, CGAL_BGL_NP_CLASS>::type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(CGAL::vertex_point, sm));

  typename Mesh::template Property_map<Vertex_index, Normal> vnormals;
  typename Mesh::template Property_map<Vertex_index, Texture> vtextures;
  typename Mesh::template Property_map<Vertex_index, Color> vcolors;
  typename Mesh::template Property_map<Face_index, Color> fcolors;
  bool created;

  std::tie(vnormals, created) = sm.template add_property_map<Vertex_index, Normal>("v:normal");
  CGAL_assertion(created);
  std::tie(vcolors, created) = sm.template add_property_map<Vertex_index, Color>("v:color", Color(0,0,0));
  CGAL_assertion(created);
  std::tie(vtextures, created) = sm.template add_property_map<Vertex_index, Texture>("v:texcoord");
  CGAL_assertion(created);
  std::tie(fcolors, created) = sm.template add_property_map<Face_index, Color>("f:color", Color(0,0,0));
  CGAL_assertion(created);

  return IO::internal::read_OFF_BGL(is, sm, CGAL::parameters::vertex_point_map(vpm)
                                                             .vertex_normal_map(vnormals)
                                                             .vertex_color_map(vcolors)
                                                             .vertex_texture_map(vtextures)
                                                             .face_color_map(fcolors));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Write

namespace IO {
namespace internal {

template <typename Point, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_with_or_without_fcolors(std::ostream& os,
                                       const Surface_mesh<Point>& sm,
                                       const CGAL_BGL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                            Mesh;
  typedef typename Mesh::Face_index                                      Face_index;
  typedef CGAL::Color                                                    Color;

  typename Mesh::template Property_map<Face_index, Color> fcolors;
  bool has_fcolors;
  std::tie(fcolors, has_fcolors) = sm.template property_map<Face_index, CGAL::Color>("f:color");

  if(has_fcolors)
    return write_OFF_BGL(os, sm, np.face_color_map(fcolors));
  else
    return write_OFF_BGL(os, sm, np);
}

template <typename Point, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_with_or_without_vtextures(std::ostream& os,
                                        const Surface_mesh<Point>& sm,
                                        const CGAL_BGL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                            Mesh;
  typedef typename Mesh::Vertex_index                                    Vertex_index;

  typedef typename Kernel_traits<Point>::Kernel                          K;
  typedef typename K::Point_2                                            Texture;

  typename Mesh::template Property_map<Vertex_index, Texture> vtextures;
  bool has_vtextures;
  std::tie(vtextures, has_vtextures) = sm.template property_map<Vertex_index, Texture>("v:texcoord");

  if(has_vtextures)
    return write_OFF_with_or_without_fcolors(os, sm, np.vertex_texture_map(vtextures));
  else
    return write_OFF_with_or_without_fcolors(os, sm, np);
}

template <typename Point, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_with_or_without_vcolors(std::ostream& os,
                                       const Surface_mesh<Point>& sm,
                                       const CGAL_BGL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                            Mesh;
  typedef typename Mesh::Vertex_index                                    Vertex_index;
  typedef CGAL::Color                                                    Color;

  typename Mesh::template Property_map<Vertex_index, Color> vcolors;
  bool has_vcolors;
  std::tie(vcolors, has_vcolors) = sm.template property_map<Vertex_index, CGAL::Color>("v:color");

  if(has_vcolors)
    return write_OFF_with_or_without_vtextures(os, sm, np.vertex_color_map(vcolors));
  else
    return write_OFF_with_or_without_vtextures(os, sm, np);
}

template <typename Point, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_with_or_without_vnormals(std::ostream& os,
                                        const Surface_mesh<Point>& sm,
                                        const CGAL_BGL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                            Mesh;
  typedef typename Mesh::Vertex_index                                    Vertex_index;

  typedef typename Kernel_traits<Point>::Kernel                          K;
  typedef typename K::Vector_3                                           Normal;

  typename Mesh::template Property_map<Vertex_index, Normal> vnormals;
  bool has_vnormals;
  std::tie(vnormals, has_vnormals) = sm.template property_map<Vertex_index, Normal>("v:normal");

  if(has_vnormals)
    return write_OFF_with_or_without_vcolors(os, sm, np.vertex_normal_map(vnormals));
  else
    return write_OFF_with_or_without_vcolors(os, sm, np);
}

} // namespace internal
} // namespace IO

/// \relates Surface_mesh
///
/// Inserts the surface mesh in an output stream.
///
/// \note The <A HREF="https://en.cppreference.com/w/cpp/io/ios_base/precision" >`precision()`</A>
///       of the output stream might not be sufficient depending on the data to be written.
///
template <typename Point, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(std::ostream& os,
               const Surface_mesh<Point>& sm,
               const CGAL_BGL_NP_CLASS& np)
{
  // Just to discard any excess named parameters
  typename CGAL::GetVertexPointMap<Surface_mesh<Point>, CGAL_BGL_NP_CLASS>::const_type
      vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                         get_const_property_map(CGAL::vertex_point, sm));

  return IO::internal::write_OFF_with_or_without_vnormals(os, sm, parameters::vertex_point_map(vpm));
}

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_IO_OFF_H
