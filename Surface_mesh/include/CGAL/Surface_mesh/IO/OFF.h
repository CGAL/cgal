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
#include <CGAL/Named_function_parameters.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Kernel_traits.h>

#include <iostream>
#include <tuple>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read
namespace IO {
namespace internal {

template<typename PolygonMesh, typename K,
         typename NamedParameters = parameters::Default_named_parameters >
class GetVertexNormalMap
{
  typedef typename PolygonMesh::template Property_map<typename PolygonMesh::Vertex_index,
                                                      typename K::Vector_3>                     DefaultMap;
  typedef DefaultMap                                                                            DefaultMap_const;

public:
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
                                                       NamedParameters, DefaultMap>::type       type;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
                                                       NamedParameters, DefaultMap_const>::type const_type;
};

template<typename PolygonMesh,
         typename NamedParameters = parameters::Default_named_parameters >
class GetVertexColorMap
{
  typedef typename PolygonMesh::template Property_map<typename PolygonMesh::Vertex_index,
                                                      CGAL::IO::Color>                          DefaultMap;
  typedef DefaultMap                                                                            DefaultMap_const;
public:
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_color_map_t,
                                                       NamedParameters, DefaultMap>::type       type;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_color_map_t,
                                                       NamedParameters, DefaultMap_const>::type const_type;
};

template<typename PolygonMesh, typename K,
         typename NamedParameters = parameters::Default_named_parameters >
class GetVertexTextureMap
{
  typedef typename PolygonMesh::template Property_map<typename PolygonMesh::Vertex_index,
                                                      typename K::Point_2>                     DefaultMap;
  typedef DefaultMap                                                                           DefaultMap_const;

public:
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_texture_map_t,
                                                       NamedParameters, DefaultMap>::type       type;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_texture_map_t,
                                                       NamedParameters, DefaultMap_const>::type const_type;
};

template<typename PolygonMesh,
         typename NamedParameters = parameters::Default_named_parameters >
class GetFaceColorMap
{
  typedef typename PolygonMesh::template Property_map<typename PolygonMesh::Face_index,
                                                      CGAL::IO::Color>                          DefaultMap;
  typedef DefaultMap                                                                            DefaultMap_const;

public:
  typedef typename internal_np::Lookup_named_param_def<internal_np::face_color_map_t,
                                                       NamedParameters, DefaultMap>::type       type;
  typedef typename internal_np::Lookup_named_param_def<internal_np::face_color_map_t,
                                                       NamedParameters, DefaultMap_const>::type const_type;
};

template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OFF_with_or_without_fcolors(std::istream& is,
                                      Surface_mesh<Point>& sm,
                                      const CGAL::File_scanner_OFF& scanner,
                                      const CGAL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                                             Mesh;
  typedef typename Mesh::Face_index                                                       Face_index;
  typedef CGAL::IO::Color                                                                 Color;

  typedef typename GetFaceColorMap<Mesh, CGAL_NP_CLASS>::type                             FCM;

  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;

  typename Mesh::template Property_map<Face_index, Color> fcm;

  bool is_fcm_requested = !(is_default_parameter<CGAL_NP_CLASS, internal_np::face_color_map_t>::value);
  if(!is_fcm_requested && scanner.has_colors())
  {
    bool created;
    std::tie(fcm, created) = sm.template add_property_map<Face_index, Color>("f:color", Color(0,0,0));
    CGAL_assertion(created);
    is_fcm_requested = true;
  }

  if(is_fcm_requested)
  {
    FCM fcolors = choose_parameter(get_parameter(np, internal_np::face_color_map), fcm);
    return CGAL::IO::internal::read_OFF_BGL(is, sm, np.face_color_map(fcolors));
  }
  else
  {
    return CGAL::IO::internal::read_OFF_BGL(is, sm, np);
  }
}

template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OFF_with_or_without_vtextures(std::istream& is,
                                        Surface_mesh<Point>& sm,
                                        const CGAL::File_scanner_OFF& scanner,
                                        const CGAL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                                                Mesh;
  typedef typename Mesh::Vertex_index                                                        Vertex_index;

  typedef typename GetK<Surface_mesh<Point>, CGAL_NP_CLASS>::Kernel                          K;
  typedef typename K::Point_2                                                                Texture;
  typedef typename GetVertexTextureMap<Mesh, K, CGAL_NP_CLASS>::type VTM;

  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;

  typename Mesh::template Property_map<Vertex_index, Texture> vtm;

  bool is_vtm_requested = !(is_default_parameter<CGAL_NP_CLASS, internal_np::vertex_texture_map_t>::value);
  if(!is_vtm_requested && scanner.has_textures())
  {
    bool created;
    std::tie(vtm, created) = sm.template add_property_map<Vertex_index, Texture>("v:texcoord");
    CGAL_assertion(created);
    is_vtm_requested = true;
  }

  if(is_vtm_requested)
  {
    VTM vtextures = choose_parameter(get_parameter(np, internal_np::vertex_texture_map), vtm);
    return read_OFF_with_or_without_fcolors(is, sm, scanner, np.vertex_texture_map(vtextures));
  }
  else
  {
    return read_OFF_with_or_without_fcolors(is, sm, scanner, np);
  }
}

template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OFF_with_or_without_vcolors(std::istream& is,
                                      Surface_mesh<Point>& sm,
                                      const CGAL::File_scanner_OFF& scanner,
                                      const CGAL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                                            Mesh;
  typedef typename Mesh::Vertex_index                                                    Vertex_index;

  typedef CGAL::IO::Color                                                                    Color;
  typedef typename GetVertexColorMap<Mesh, CGAL_NP_CLASS>::type  VCM;

  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;

  typename Mesh::template Property_map<Vertex_index, Color> vcm;

  bool is_vcm_requested = !(is_default_parameter<CGAL_NP_CLASS,  internal_np::vertex_color_map_t>::value);
  if(!is_vcm_requested && scanner.has_colors())
  {
    bool created;
    std::tie(vcm, created) = sm.template add_property_map<Vertex_index, Color>("v:color", Color(0,0,0));
    CGAL_assertion(created);
    is_vcm_requested = true;
  }

  if(is_vcm_requested)
  {
    VCM vcolors = choose_parameter(get_parameter(np, internal_np::vertex_color_map), vcm);
    return read_OFF_with_or_without_vtextures(is, sm, scanner, np.vertex_color_map(vcolors));
  }
  else
  {
    return read_OFF_with_or_without_vtextures(is, sm, scanner, np);
  }
}

template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OFF_with_or_without_vnormals(std::istream& is,
                                       Surface_mesh<Point>& sm,
                                       const CGAL::File_scanner_OFF& scanner,
                                       const CGAL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                                                Mesh;
  typedef typename Mesh::Vertex_index                                                        Vertex_index;

  typedef typename GetK<Surface_mesh<Point>, CGAL_NP_CLASS>::Kernel                      K;
  typedef typename K::Vector_3                                                               Normal;
  typedef typename GetVertexNormalMap<Mesh, K, CGAL_NP_CLASS>::type  VNM;

  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;

  typename Mesh::template Property_map<Vertex_index, Normal> vnm;

  bool is_vnm_requested = !(is_default_parameter<CGAL_NP_CLASS, internal_np::vertex_normal_map_t>::value);
  if(!is_vnm_requested && scanner.has_normals())
  {
    bool created;
    std::tie(vnm, created) = sm.template add_property_map<Vertex_index, Normal>("v:normal");
    CGAL_assertion(created);
    is_vnm_requested = true;
  }

  if(is_vnm_requested)
  {
    VNM vnormals = choose_parameter(get_parameter(np, internal_np::vertex_normal_map), vnm);
    return read_OFF_with_or_without_vcolors(is, sm, scanner, np.vertex_normal_map(vnormals));
  }
  else
  {
    return read_OFF_with_or_without_vcolors(is, sm, scanner, np);
  }
}

} // namespace internal

/// \ingroup PkgSurfaceMeshIOFuncOFF
///
/// \brief extracts the surface mesh from an input stream in the \ref IOStreamOFF
///        and appends it to the surface mesh `sm`.
///
/// This function reads points, as well as vertex normals, vertex and face colors,
/// and texture vertex coordinates if those attributes are available in the input.
/// These last four attributes are stored in internal property maps of `sm`
/// named "v:normal", "v:color", "f:color", and `"v:texcoord"`, respectively,
/// which will be created if they do not already exist.
/// If property maps are passed through named parameters (see below),
/// then they are used instead of the internal ones.
///
/// Ignores comment lines which start with a hash, and lines with whitespace.
///
/// \attention The graph `sm` is not cleared, and the data from the stream is added.
///
/// \tparam Point The type of the \em point property of a vertex. There is no requirement on `P`,
///               besides being default constructible and assignable.
///               In typical use cases it will be a 2D or 3D point type.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param is the input stream
/// \param sm the surface mesh to be constructed
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `sm`}
///     \cgalParamType{a class model of `WritablePropertyMap` with `Surface_mesh::Vertex_index`
///                    as key type and `Point` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{vertex_normal_map}
///     \cgalParamDescription{a property map associating normals to the vertices of `sm`}
///     \cgalParamType{a class model of `WritablePropertyMap` with `Surface_mesh::Vertex_index`
///                    as key type and a 3D vector type issued from the same kernel as `Point` as value type}
///     \cgalParamDefault{If this parameter is unused, vertex normals (if they exist)
///                       will be written in an internal property map called `v:normal`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{vertex_color_map}
///     \cgalParamDescription{a property map associating colors to the vertices of `sm`}
///     \cgalParamType{a class model of `WritablePropertyMap` with `Surface_mesh::Vertex_index`
///                    as key type and `CGAL::IO::Color` as value type}
///     \cgalParamDefault{If this parameter is unused, vertex colors (if they exist)
///                       will be written in an internal property map called `v:color`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{vertex_texture_map}
///     \cgalParamDescription{a property map associating textures to the vertices of `sm`}
///     \cgalParamType{a class model of `WritablePropertyMap` with `Surface_mesh::Vertex_index`
///                    as key type and a 2D vector type issued from the same kernel as `Point` as value type}
///     \cgalParamDefault{If this parameter is unused, vertex textures (if they exist)
///                       will be written in an internal property map called `v:texcoords`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{face_color_map}
///     \cgalParamDescription{a property map associating colors to the faces of `sm`}
///     \cgalParamType{a class model of `WritablePropertyMap` with `Surface_mesh::Face_index`
///                    as key type and `CGAL::IO::Color` as value type}
///     \cgalParamDefault{If this parameter is unused, face colors (if they exist)
///                       will be written in an internal property map called `f:color`.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \pre The data in the stream must represent a two-manifold. If this is not the case
///      the `failbit` of `is` is set and the mesh cleared.
///
template <typename Point,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(std::istream& is,
              Surface_mesh<Point>& sm,
              const CGAL_NP_CLASS& np = parameters::default_values())
{
  std::streampos pos = is.tellg();
  CGAL::File_scanner_OFF scanner(is, false);
  is.seekg(pos);

  bool res = internal::read_OFF_with_or_without_vnormals(is, sm, scanner, np);
  if(!res)
    sm.clear();

  return res;
}

template <typename Point,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const std::string& fname,
              Surface_mesh<Point>& sm,
              const CGAL_NP_CLASS& np = parameters::default_values())
{
  std::ifstream in(fname.c_str());
  return read_OFF(in, sm, np);
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgSurfaceMeshIOFuncDeprecated
  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::read_OFF(std::istream&, const Surface_mesh<Point>&)` should be used instead.
*/
template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool read_off(std::istream& is, Surface_mesh<Point>& sm, const CGAL_NP_CLASS& np = parameters::default_values())
{
  return IO::read_OFF(is, sm, np);
}

/*!
  \ingroup PkgSurfaceMeshIOFuncDeprecated
  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::read_OFF(std::istream&, const Surface_mesh<Point>&)` should be used instead.
*/
template <typename Point>
CGAL_DEPRECATED bool read_off(Surface_mesh<Point>& sm, const std::string& filename)
{
  return IO::read_OFF(filename, sm, parameters::default_values());
}
#endif // CGAL_NO_DEPRECATED_CODE

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace IO {
namespace internal {

template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_with_or_without_fcolors(std::ostream& os,
                                       const Surface_mesh<Point>& sm,
                                       const CGAL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                            Mesh;
  typedef typename Mesh::Face_index                                      Face_index;

  using parameters::is_default_parameter;

  const bool has_fcolors = !(is_default_parameter<CGAL_NP_CLASS, internal_np::face_color_map_t>::value);

  auto fcolors  = sm.template property_map<Face_index, CGAL::IO::Color>("f:color");

  if(!has_fcolors && fcolors.has_value() && std::distance(fcolors.value().begin(), fcolors.value().end()) > 0)
    return write_OFF_BGL(os, sm, np.face_color_map(fcolors.value()));
  else
    return write_OFF_BGL(os, sm, np);
}

template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_with_or_without_vtextures(std::ostream& os,
                                         const Surface_mesh<Point>& sm,
                                         const CGAL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                            Mesh;
  typedef typename Mesh::Vertex_index                                    Vertex_index;

  typedef typename GetK<Surface_mesh<Point>, CGAL_NP_CLASS>::Kernel      K;
  typedef typename K::Point_2                                            Texture;

  using parameters::is_default_parameter;

  const bool has_vtextures = !(is_default_parameter<CGAL_NP_CLASS, internal_np::vertex_texture_map_t>::value);

  auto vtextures = sm.template property_map<Vertex_index, Texture>("v:texcoord");

  if(!has_vtextures && vtextures.has_value() && std::distance(vtextures.value().begin(), vtextures.value().end()) > 0)
    return write_OFF_with_or_without_fcolors(os, sm, np.vertex_texture_map(vtextures.value()));
  else
    return write_OFF_with_or_without_fcolors(os, sm, np);
}

template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_with_or_without_vcolors(std::ostream& os,
                                       const Surface_mesh<Point>& sm,
                                       const CGAL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                            Mesh;
  typedef typename Mesh::Vertex_index                                    Vertex_index;

  using parameters::is_default_parameter;

  const bool has_vcolors = !(is_default_parameter<CGAL_NP_CLASS, internal_np::vertex_color_map_t>::value);


  auto vcolors = sm.template property_map<Vertex_index, CGAL::IO::Color>("v:color");

  if(!has_vcolors && vcolors.has_value() && std::distance(vcolors.value().begin(), vcolors.value().end()) > 0)
    return write_OFF_with_or_without_vtextures(os, sm, np.vertex_color_map(vcolors.value()));
  else
    return write_OFF_with_or_without_vtextures(os, sm, np);
}

template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OFF_with_or_without_vnormals(std::ostream& os,
                                        const Surface_mesh<Point>& sm,
                                        const CGAL_NP_CLASS& np)
{
  typedef Surface_mesh<Point>                                            Mesh;
  typedef typename Mesh::Vertex_index                                    Vertex_index;

  typedef typename GetK<Surface_mesh<Point>, CGAL_NP_CLASS>::Kernel      K;
  typedef typename K::Vector_3                                           Normal;

  using parameters::is_default_parameter;

  const bool has_vnormals = !(is_default_parameter<CGAL_NP_CLASS, internal_np::vertex_normal_map_t>::value);

  auto vnormals = sm.template property_map<Vertex_index, Normal>("v:normal");

  if(!has_vnormals && vnormals.has_value() && std::distance(vnormals.value().begin(), vnormals.value().end()) > 0)
    return write_OFF_with_or_without_vcolors(os, sm, np.vertex_normal_map(vnormals.value()));
  else
    return write_OFF_with_or_without_vcolors(os, sm, np);
}

} // namespace internal

/// \ingroup PkgSurfaceMeshIOFuncOFF
///
/// \brief writes the surface mesh `sm` in the output stream, using the \ref IOStreamOFF.
///
/// This overload of \link PkgBGLIOFct `write_OFF(std::ostream&, const Graph&)` \endlink will also output
/// the following property maps internal to the surface mesh, if they exist and if they are not
/// already present in the named parameters:
///
/// - vertex normals (property map named "v:normal" in the surface mesh)
/// - vertex colors (property map named "v:color" in the surface mesh)
/// - vertex textures (property map named "v:texcoord" in the surface mesh)
/// - face colors (property map named "f:color" in the surface mesh)
///
/// \tparam Point The type of the \em point property of a vertex. There is no requirement on `P`,
///               besides being default constructible and assignable.
///               In typical use cases it will be a 2D or 3D point type.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param os the output stream
/// \param sm the surface mesh to be output
/// \param np optional \ref bgl_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `sm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `Surface_mesh::Vertex_index`
///                    as key type and `%Point` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{vertex_normal_map}
///     \cgalParamDescription{a property map associating normals to the vertices of `sm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `Surface_mesh::Vertex_index`
///                    as key type and a 3D vector type issued from the same kernel as `Point` as value type}
///     \cgalParamDefault{vertex normals will be output using the internal property map, if it exists.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{vertex_color_map}
///     \cgalParamDescription{a property map associating colors to the vertices of `sm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `Surface_mesh::Vertex_index`
///                    as key type and `CGAL::IO::Color` as value type}
///     \cgalParamDefault{vertex colors will be output using the internal property map, if it exists.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{vertex_texture_map}
///     \cgalParamDescription{a property map associating textures to the vertices of `sm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `Surface_mesh::Vertex_index`
///                    as key type and a 2D point type issued from the same kernel as `Point` as value type}
///     \cgalParamDefault{vertex textures will be output using the internal property map, if it exists.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{face_color_map}
///     \cgalParamDescription{a property map associating colors to the faces of `sm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `Surface_mesh::Face_index`
///                    as key type and `CGAL::IO::Color` as value type}
///     \cgalParamDefault{face colors will be output using the internal property map, if it exists.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{stream_precision}
///     \cgalParamDescription{a parameter used to set the precision (i.e. how many digits are generated) of the output stream}
///     \cgalParamType{int}
///     \cgalParamDefault{the precision of the stream `os`}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \returns `true` if writing was successful, `false` otherwise.
///
template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(std::ostream& os,
               const Surface_mesh<Point>& sm,
               const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::is_default_parameter;

  const bool has_vpoints = !(is_default_parameter<CGAL_NP_CLASS, internal_np::vertex_point_t>::value);
  if(has_vpoints)
    return internal::write_OFF_with_or_without_vnormals(os, sm, np);

  return internal::write_OFF_with_or_without_vnormals(os, sm, np.vertex_point_map(get_const_property_map(CGAL::vertex_point, sm)));
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE

/*!
  \ingroup PkgSurfaceMeshIOFuncDeprecated
  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::write_OFF(std::ostream&, const Surface_mesh<Point>&)` should be used instead.
*/
template <typename Point, typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL_DEPRECATED bool write_off(std::ostream& os, const Surface_mesh<Point>& sm, const CGAL_NP_CLASS& np = parameters::default_values())
{
  return IO::write_OFF(os, sm, np);
}

/*!
  \ingroup PkgSurfaceMeshIOFuncDeprecated
  \deprecated This function is deprecated since \cgal 5.3, `CGAL::IO::write_OFF(std::ostream&, const Surface_mesh<Point>&)` should be used instead.
*/
template <typename Point>
CGAL_DEPRECATED bool write_off(const Surface_mesh<Point>& sm, const std::string& filename)
{
  return IO::write_OFF(filename, sm, parameters::default_values());
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_IO_OFF_H
