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

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#endif

namespace CGAL {
namespace internal{
namespace IO{
template<typename PolygonMesh,
         typename K,
         typename NamedParameters = Named_function_parameters<bool, internal_np::all_default_t> >
class GetVertexNormalMap
{
  typedef typename PolygonMesh::template Property_map<typename PolygonMesh::Vertex_index, typename K::Vector_3>
  DefaultMap;
  typedef DefaultMap DefaultMap_const;
public:
  typedef typename internal_np::Lookup_named_param_def<
  internal_np::vertex_normal_map_t,
  NamedParameters,
  DefaultMap
  > ::type  type;
  typedef typename internal_np::Lookup_named_param_def<
    internal_np::vertex_normal_map_t,
    NamedParameters,
    DefaultMap_const
    > ::type  const_type;
};

template<typename PolygonMesh,
         typename NamedParameters = Named_function_parameters<bool, internal_np::all_default_t> >
class GetVertexColorMap
{
  typedef typename PolygonMesh::template Property_map<typename PolygonMesh::Vertex_index, CGAL::Color>
  DefaultMap;
  typedef DefaultMap DefaultMap_const;
public:
  typedef typename internal_np::Lookup_named_param_def<
  internal_np::vertex_color_map_t,
  NamedParameters,
  DefaultMap
  > ::type  type;
  typedef typename internal_np::Lookup_named_param_def<
    internal_np::vertex_color_map_t,
    NamedParameters,
    DefaultMap_const
    > ::type  const_type;
};

template<typename PolygonMesh,
         typename K,
         typename NamedParameters = Named_function_parameters<bool, internal_np::all_default_t> >
class GetVertexTextureMap
{
  typedef typename PolygonMesh::template Property_map<typename PolygonMesh::Vertex_index, typename K::Point_2>
  DefaultMap;
  typedef DefaultMap DefaultMap_const;
public:
  typedef typename internal_np::Lookup_named_param_def<
  internal_np::vertex_texture_map_t,
  NamedParameters,
  DefaultMap
  > ::type  type;
  typedef typename internal_np::Lookup_named_param_def<
    internal_np::vertex_texture_map_t,
    NamedParameters,
    DefaultMap_const
    > ::type  const_type;
};
template<typename PolygonMesh,
         typename NamedParameters = Named_function_parameters<bool, internal_np::all_default_t> >
class GetFaceColorMap
{
  typedef typename PolygonMesh::template Property_map<typename PolygonMesh::Face_index, CGAL::Color>
  DefaultMap;
  typedef DefaultMap DefaultMap_const;
public:
  typedef typename internal_np::Lookup_named_param_def<
  internal_np::face_color_map_t,
  NamedParameters,
  DefaultMap
  > ::type  type;
  typedef typename internal_np::Lookup_named_param_def<
    internal_np::face_color_map_t,
    NamedParameters,
    DefaultMap_const
    > ::type  const_type;
};
}//end IO
}// end internal
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Read

/// \relates Surface_mesh
/// \ingroup PkgSurfaceMeshIOFunc
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
template <typename Point,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(std::istream& is,
              Surface_mesh<Point>& sm,
              const CGAL_BGL_NP_CLASS& np,
              bool verbose = true)
{
  typedef Surface_mesh<Point>                                            Mesh;
  typedef typename Mesh::Vertex_index                                    Vertex_index;
  typedef typename Mesh::Face_index                                      Face_index;

  typedef typename GetK<Surface_mesh<Point>, CGAL_BGL_NP_CLASS>::Kernel  K;
  typedef typename K::Vector_3                                           Normal;
  typedef typename K::Point_2                                            Texture;
  typedef CGAL::Color                                                    Color;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;


  typedef typename CGAL::GetVertexPointMap<Mesh, CGAL_BGL_NP_CLASS>::type   VPM;
  typedef typename
  CGAL::internal::IO::GetVertexNormalMap<Mesh, K, CGAL_BGL_NP_CLASS>::type  VNM;
typedef typename
  CGAL::internal::IO::GetVertexColorMap<Mesh, CGAL_BGL_NP_CLASS>::type      VCM;
  typedef typename
  CGAL::internal::IO::GetVertexTextureMap<Mesh, K, CGAL_BGL_NP_CLASS>::type VTM;
  typedef typename
  CGAL::internal::IO::GetFaceColorMap<Mesh, CGAL_BGL_NP_CLASS>::type        FCM;

  const bool is_vnm_requested =
      !(is_default_parameter(get_parameter(np, internal_np::vertex_normal_map)));
  const bool is_vcm_requested =
      !(is_default_parameter(get_parameter(np, internal_np::vertex_color_map)));
  const bool is_vtm_requested =
      !(is_default_parameter(get_parameter(np, internal_np::vertex_texture_map)));
  const bool is_fcm_requested =
      !(is_default_parameter(get_parameter(np, internal_np::face_color_map)));

  bool created;
  typename Mesh::template Property_map<Vertex_index, Normal> vnm;
  typename Mesh::template Property_map<Vertex_index, Color> vcm;
  typename Mesh::template Property_map<Vertex_index, Texture> vtm;
  typename Mesh::template Property_map<Face_index, Color> fcm;

  if(!is_vnm_requested)
  {
    std::tie(vnm, created) =
        sm.template add_property_map<Vertex_index, Normal>("v:normal");
    CGAL_assertion(created);
  }
  if(!is_vcm_requested)
  {
    std::tie(vcm, created) =
        sm.template add_property_map<Vertex_index, Color>("v:color", Color(0,0,0));
    CGAL_assertion(created);
  }
  if(!is_vtm_requested)
  {
    std::tie(vtm, created) =
        sm.template add_property_map<Vertex_index, Texture>("v:texcoord");
    CGAL_assertion(created);
  }
  if(!is_fcm_requested)
  {
    std::tie(fcm, created) =
        sm.template add_property_map<Face_index, Color>("f:color", Color(0,0,0));
    CGAL_assertion(created);
  }

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(CGAL::vertex_point, sm));
  VNM vnormals = choose_parameter(get_parameter(np, internal_np::vertex_normal_map),
                                  vnm);
  VCM vcolors = choose_parameter(get_parameter(np, internal_np::vertex_color_map),
                                 vcm);
  VTM vtextures = choose_parameter(get_parameter(np, internal_np::vertex_texture_map),
                                   vtm);
  FCM fcolors = choose_parameter(get_parameter(np, internal_np::face_color_map),
                                 fcm);

  bool res =  IO::internal::read_OFF_BGL(is, sm, CGAL::parameters::vertex_point_map(vpm)
                                         .vertex_normal_map(vnormals)
                                         .vertex_color_map(vcolors)
                                         .vertex_texture_map(vtextures)
                                         .face_color_map(fcolors),
                                         verbose);
  if(!res)
    sm.clear();
  return res;
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
  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;
  const bool has_fcolors = !(is_default_parameter(get_parameter(np, internal_np::face_color_map)));
  bool has_internal_fcolors;
  std::tie(fcolors, has_internal_fcolors) = sm.template property_map<Face_index, CGAL::Color>("f:color");

  if(!has_fcolors && has_internal_fcolors && !std::distance(fcolors.begin(), fcolors.end()) > 0)
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

  typedef typename GetK<Surface_mesh<Point>, CGAL_BGL_NP_CLASS>::Kernel  K;
  typedef typename K::Point_2                                            Texture;

  typename Mesh::template Property_map<Vertex_index, Texture> vtextures;
  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;
  const bool has_vtextures = !(is_default_parameter(get_parameter(np, internal_np::vertex_texture_map)));
  bool has_internal_vtextures;
  std::tie(vtextures, has_internal_vtextures) = sm.template property_map<Vertex_index, Texture>("v:texcoord");

  if(!has_vtextures && has_internal_vtextures && !std::distance(vtextures.begin(), vtextures.end()) > 0)
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
  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;
  const bool has_vcolors = !(is_default_parameter(get_parameter(np, internal_np::vertex_color_map)));
  bool has_internal_vcolors;


  std::tie(vcolors, has_internal_vcolors) = sm.template property_map<Vertex_index, CGAL::Color>("v:color");

  if(!has_vcolors && has_internal_vcolors && !std::distance(vcolors.begin(), vcolors.end()) > 0)
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

  typedef typename GetK<Surface_mesh<Point>, CGAL_BGL_NP_CLASS>::Kernel  K;
  typedef typename K::Vector_3                                           Normal;

  typename Mesh::template Property_map<Vertex_index, Normal> vnormals;
  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;
  const bool has_vnormals = !(is_default_parameter(get_parameter(np, internal_np::vertex_normal_map)));
  bool has_internal_vnormals;
  std::tie(vnormals, has_internal_vnormals) = sm.template property_map<Vertex_index, Normal>("v:normal");

  if(!has_vnormals && has_internal_vnormals && !std::distance(vnormals.begin(), vnormals.end()) > 0)
    return write_OFF_with_or_without_vcolors(os, sm, np.vertex_normal_map(vnormals));
  else
    return write_OFF_with_or_without_vcolors(os, sm, np);
}

} // namespace internal
} // namespace IO

/// \relates Surface_mesh
/// \ingroup PkgSurfaceMeshIOFunc
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
  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;

  const bool has_vpoints= !(is_default_parameter(get_parameter(np, internal_np::vertex_point)));
  if(has_vpoints)
    return IO::internal::write_OFF_with_or_without_vnormals(os, sm, np);
  return IO::internal::write_OFF_with_or_without_vnormals(os, sm, np.vertex_point_map(get_const_property_map(CGAL::vertex_point, sm)));
}

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_IO_OFF_H
