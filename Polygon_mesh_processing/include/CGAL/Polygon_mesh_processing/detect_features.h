// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Stephane Tayeb, Maxime Gimeno
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_H
#define CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_H

#include <CGAL/license/Polygon_mesh_processing/detect_features.h>

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <set>

namespace CGAL {
namespace Polygon_mesh_processing {

namespace internal
{
template <typename Int>
Int
generate_patch_id(Int, int i)
{
  return Int(i);
}

template <typename Int>
std::pair<Int, Int>
generate_patch_id(std::pair<Int, Int>, int i)
{
  return std::pair<Int, Int>(i, 0);
}

template <typename PolygonMesh, typename VPM, typename GT>
bool
is_sharp(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
         const PolygonMesh& pmesh,
         const VPM vpm,
         const GT gt,
         const typename GT::Sign cos_sign,
         const typename GT::FT sq_cos_angle)
{
  typedef typename GT::FT                                                      FT;
  typedef typename GT::Vector_3                                                Vector_3;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;

  if(is_border_edge(h, pmesh))
    return false;

  typename GT::Construct_vector_3 vector = gt.construct_vector_3_object();
  typename GT::Construct_cross_product_vector_3 cross = gt.construct_cross_product_vector_3_object();
  typename GT::Compute_squared_length_3 sq_length = gt.compute_squared_length_3_object();

  halfedge_descriptor opp_h = opposite(h, pmesh);

  const Vector_3 vc = vector(get(vpm, source(h, pmesh)), get(vpm, target(h, pmesh)));
  const Vector_3 v1 = vector(get(vpm, source(h, pmesh)), get(vpm, target(next(h, pmesh), pmesh)));
  const Vector_3 v2 = vector(get(vpm, target(next(opp_h, pmesh), pmesh)), get(vpm, target(opp_h, pmesh)));

  const Vector_3 n1 = cross(vc, v1);
  const Vector_3 n2 = cross(vc, v2);
  const FT sp = n1 * n2;

  // n1.n2 <= cos() without computing the norms of the vectors
  if(cos_sign == NEGATIVE)
    return (is_negative(sp) && (square(sp) >= sq_cos_angle * sq_length(n1) * sq_length(n2)));
  else
    return (is_negative(sp) || (square(sp) <= sq_cos_angle * sq_length(n1) * sq_length(n2)));
}

//wrapper for patchid map.
template<typename PatchIdMap,
         typename ValueType = typename boost::property_traits<PatchIdMap>::value_type>
struct PatchIdMapWrapper
{
  typedef typename boost::read_write_property_map_tag category;
  typedef ValueType value_type;
  typedef typename boost::property_traits<PatchIdMap>::reference reference;
  typedef typename boost::property_traits<PatchIdMap>::key_type key_type;

  PatchIdMap map;
  int offset;
  PatchIdMapWrapper(PatchIdMap map, int offset)
    : map(map), offset(offset){}
};

template <typename PatchIdMap, typename Handle_type, typename Int>
typename PatchIdMapWrapper<PatchIdMap, Int>::value_type
get(PatchIdMapWrapper<PatchIdMap, Int>& map, Handle_type h)
{
  typedef typename PatchIdMapWrapper<PatchIdMap, Int>::value_type value_type;
  return value_type(get(map.map, h) - map.offset);
}

template <typename PatchIdMap, typename Handle_type, typename Int>
void put(PatchIdMapWrapper<PatchIdMap, Int>& map, Handle_type h,
         typename PatchIdMapWrapper<PatchIdMap, Int>::value_type pid)
{
  typedef typename PatchIdMapWrapper<PatchIdMap, Int>::value_type value_type;
  put(map.map, h, value_type(pid + map.offset));
}


//specialization for std::pair
template<typename PatchIdMap, typename Int>
struct PatchIdMapWrapper<PatchIdMap, std::pair<Int, Int> >
{
  typedef typename boost::read_write_property_map_tag category;
  typedef Int value_type;
  typedef typename boost::property_traits<PatchIdMap>::reference reference;
  typedef typename boost::property_traits<PatchIdMap>::key_type key_type;

  PatchIdMap map;
  int offset;
  PatchIdMapWrapper(PatchIdMap map, int offset)
    : map(map), offset(offset){}
};

template <typename PatchIdMap, typename Handle_type, typename Int>
typename PatchIdMapWrapper<PatchIdMap, std::pair<Int, Int> >::value_type
get(PatchIdMapWrapper<PatchIdMap, std::pair<Int, Int> >& map, Handle_type h)
{
  return Int(get(map.map, h).first - map.offset);
}

template <typename PatchIdMap, typename Handle_type, typename Int>
void put(PatchIdMapWrapper<PatchIdMap, std::pair<Int, Int> >& map, Handle_type h,
         typename PatchIdMapWrapper<PatchIdMap, std::pair<Int, Int> >::value_type pid)
{
  put(map.map, h, std::pair<Int, Int>(Int(pid+map.offset), 0));
}

template <typename PolygonMesh, typename PatchIdMap,
          typename EdgeIsFeatureMap, typename NamedParameters = parameters::Default_named_parameters>
typename boost::graph_traits<PolygonMesh>::faces_size_type
detect_surface_patches(const PolygonMesh& p,
                       PatchIdMap patch_id_map,
                       EdgeIsFeatureMap eif,
                       const NamedParameters& np = parameters::default_values())
{
  int offset = static_cast<int>(
          parameters::choose_parameter(parameters::get_parameter(np, internal_np::first_index), 1));

  internal::PatchIdMapWrapper<PatchIdMap,
                              typename boost::property_traits<PatchIdMap>::value_type>
          wrapmap(patch_id_map, offset);

  return connected_components(p, wrapmap,
                              parameters::edge_is_constrained_map(eif)
                                         .face_index_map(CGAL::get_initialized_face_index_map(p, np)));
}

template <typename FT, typename PolygonMesh, typename VPM, typename GT, typename EIFMap, typename VNFEMap>
void sharp_call(const FT angle_in_deg,
                const PolygonMesh& pmesh,
                const VPM vpm,
                const GT gt,
                EIFMap edge_is_feature_map,
                VNFEMap vnfe)
{
  // Initialize vertices
  for(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd : vertices(pmesh))
    put(vnfe, vd, 0);

  const FT cos_angle = std::cos(CGAL::to_double(angle_in_deg) * CGAL_PI / 180.);
  const FT sq_cos_angle = square(cos_angle);

  // Detect sharp edges
  for(typename boost::graph_traits<PolygonMesh>::edge_descriptor ed : edges(pmesh))
  {
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor he = halfedge(ed, pmesh);
    if(is_border_edge(he, pmesh) ||
       angle_in_deg == FT() ||
       (angle_in_deg != FT(180) && internal::is_sharp(he, pmesh, vpm, gt, CGAL::sign(cos_angle), sq_cos_angle)))
    {
      put(edge_is_feature_map, edge(he, pmesh), true);
      put(vnfe, target(he, pmesh), get(vnfe, target(he, pmesh))+1);
      put(vnfe, source(he, pmesh), get(vnfe, source(he, pmesh))+1);
    }
  }
}

template <typename FT, typename PolygonMesh, typename VPM, typename GT, typename EIFMap>
void sharp_call(const FT angle_in_deg,
                const PolygonMesh& pmesh,
                const VPM vpm,
                const GT gt,
                EIFMap edge_is_feature_map,
                const internal_np::Param_not_found&)
{
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor     edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  const FT cos_angle = std::cos(CGAL::to_double(angle_in_deg) * CGAL_PI / 180.);
  const FT sq_cos_angle = square(cos_angle);

  // Detect sharp edges
  for(edge_descriptor ed : edges(pmesh))
  {
    halfedge_descriptor he = halfedge(ed, pmesh);
    if(is_border_edge(he, pmesh) ||
       angle_in_deg == FT() ||
       (angle_in_deg != FT(180) && internal::is_sharp(he, pmesh, vpm, gt, CGAL::sign(cos_angle), sq_cos_angle)))
    {
      put(edge_is_feature_map, edge(he, pmesh), true);
    }
  }
}

} // namespace internal

/*!
 * \ingroup PMP_detect_features_grp
 *
 * \brief detects and marks the edges that are considered to be sharp with respect to the given angle bound.
 *
 * `angle_in_deg` gives the maximum angle (in degrees) between the two normal vectors of adjacent triangles.
 * For an edge of the input polygon mesh, if the angle between the two normal vectors of its incident facets is bigger
 * than the given bound, then the edge is marked as being a feature edge.
 *
 * Also computes the number of sharp edges incident to each vertex, if `vertex_feature_degree_map` is provided.
 *
 * \tparam PolygonMesh a model of `HalfedgeListGraph`
 * \tparam FT a number type. It is
 * either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
 * or from the geometric traits class deduced from the point property map
 * of `PolygonMesh`.
 * \tparam EdgeIsFeatureMap a model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
 *  as key type and `bool` as value type. It must be default constructible.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param pmesh the polygon mesh
 * \param angle_in_deg the dihedral angle bound
 * \param edge_is_feature_map the property map that will contain the sharp-or-not status of each edge of `pmesh`
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_feature_degree_map}
 *     \cgalParamDescription{a property map that will associate to each vertex of `pmesh` the number of incident feature edges}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and `int` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_feature_degree_t(), pmesh)`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \see `sharp_edges_segmentation()`
 */
#ifdef DOXYGEN_RUNNING
template <typename PolygonMesh, typename FT,
          typename EdgeIsFeatureMap, typename NamedParameters>
#else
template <typename PolygonMesh, typename EdgeIsFeatureMap, typename NamedParameters = parameters::Default_named_parameters>
#endif
void detect_sharp_edges(const PolygonMesh& pmesh,
#ifdef DOXYGEN_RUNNING
                        FT angle_in_deg,
#else
                        typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT angle_in_deg,
#endif
                        EdgeIsFeatureMap edge_is_feature_map,
                        const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // extract types from NPs
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(boost::vertex_point, pmesh));

  internal::sharp_call(angle_in_deg, pmesh, vpm, gt, edge_is_feature_map,
                       get_parameter(np, internal_np::vertex_feature_degree));
}

/*!
 * \ingroup PMP_detect_features_grp
 *
 * collects the surface patches of the faces incident to each vertex of the input polygon mesh.
 *
 * \tparam PolygonMesh a model of `HalfedgeListGraph`
 * \tparam PatchIdMap a model of `ReadablePropertyMap` with
 *  `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type
 *  and the desired patch id, model of `CopyConstructible` as value type.
 * \tparam VertexIncidentPatchesMap a model of mutable `LvaluePropertyMap` with
 *  `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type. Its value type
 *  must be a container of `boost::property_traits<PatchIdMap>::%value_type` and have a function `insert()`.
 *  A `std::set` or a `boost::unordered_set` are recommended, as a patch index may be
 *  inserted several times.
 * \tparam EdgeIsFeatureMap a model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
 *  as key type and `bool` as value type.
 *
 * \param pmesh the polygon mesh
 * \param patch_id_map the property map containing the surface patch ids for the faces of `pmesh`. It must be already filled.
 * \param vertex_incident_patches_map a property map that will contain the patch ids of all the faces incident to each vertex of `pmesh`.
 * \param edge_is_feature_map a filled property map that will contain the sharp-or-not status of each edge of `pmesh`
 *
 * \see `sharp_edges_segmentation()`
 */

template <typename PolygonMesh, typename PatchIdMap,
          typename VertexIncidentPatchesMap, typename EdgeIsFeatureMap>
void detect_vertex_incident_patches(const PolygonMesh& pmesh,
                             const PatchIdMap patch_id_map,
                             VertexIncidentPatchesMap vertex_incident_patches_map,
                             const EdgeIsFeatureMap edge_is_feature_map)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  for(vertex_descriptor vit :vertices(pmesh))
  {
    // Look only at feature vertices
    if(!get(edge_is_feature_map, edge(halfedge(vit, pmesh), pmesh)))
      continue;

    // Loop on incident facets of vit
    typename VertexIncidentPatchesMap::value_type& id_set = vertex_incident_patches_map[vit];
    for(halfedge_descriptor he : halfedges_around_target(vit, pmesh))
    {
      if(!is_border(he, pmesh))
      {
        id_set.insert(get(patch_id_map, face(he, pmesh)));
      }
      else if(!is_border(opposite(he, pmesh), pmesh))
      {
        id_set.insert(get(patch_id_map, face(opposite(he, pmesh), pmesh)));
      }
    }
  }
}

namespace internal {

template<typename PolygonMesh, typename PIDMap, typename VIPMap, typename EIFMap>
void vip_call(const PolygonMesh& mesh, PIDMap pid, VIPMap vip, EIFMap eif)
{
  CGAL::Polygon_mesh_processing::detect_vertex_incident_patches(mesh, pid, vip, eif);
}

template<typename PolygonMesh, typename PIDMap, typename EIFMap>
void vip_call(const PolygonMesh&, PIDMap, const internal_np::Param_not_found&, EIFMap)
{
  //do nothing when the parameter is not given
}

} // namespace internal

/*!
 * \ingroup PMP_detect_features_grp
 *
 * This function calls successively `CGAL::Polygon_mesh_processing::detect_sharp_edges()`,
 * `CGAL::Polygon_mesh_processing::connected_components()`, and
 * `CGAL::Polygon_mesh_processing::detect_vertex_incident_patches()`
 *
 * It detects and marks the sharp edges of `pmesh` according to `angle_in_deg`.
 * The sharp edges define a segmentation of `pmesh`, that is done by computing a
 * surface patch id for each face.
 *
 * \tparam PolygonMesh a model of `FaceGraph`
 * \tparam FT a number type. It is
 * either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
 * or from the geometric traits class deduced from the point property map
 * of `PolygonMesh`.
 * \tparam EdgeIsFeatureMap a model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
 * \tparam PatchIdMap a model of `ReadWritePropertyMap` with
 *  `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type
 *  and the desired patch id, model of `CopyConstructible` as value type.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param pmesh the polygon mesh
 * \param angle_in_deg the dihedral angle bound
 * \param edge_is_feature_map the property map that will contain the sharp-or-not status of each edge of `pmesh`
 * \param patch_id_map the property map that will contain the surface patch ids for the faces of `pmesh`.
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_feature_degree_map}
 *     \cgalParamDescription{a property map that will associate to each vertex of `pmesh` the number of incident feature edges}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `int` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_feature_degree_t(), pmesh)`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{first_index}
 *     \cgalParamDescription{the index of the first surface patch of `pmesh`}
 *     \cgalParamType{`std::size_t`}
 *     \cgalParamExtra{The patches will be numbered on `[first_index; first_index + num_patches]`,
 *                     where `num_patches` is the number of surface patches.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{face_index_map}
 *     \cgalParamDescription{a property map associating to each face of `pmesh` a unique index between `0` and `num_faces(pmesh) - 1`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{an automatically indexed internal map}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_incident_patches_map}
 *     \cgalParamDescription{a property map that will contain the patch ids of all the faces incident to each vertex of `pmesh`}
 *     \cgalParamType{a model of mutable `LvaluePropertyMap` with
 *                    `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type. Its value type
 *                    must be a container of `boost::property_traits<PatchIdMap>::%value_type` and have a function `insert()`.}
 *     \cgalParamExtra{A `std::set` or a `boost::unordered_set` are recommended, as a patch index may be inserted several times.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 * \returns the number of surface patches.
 *
 * @see `CGAL::Polygon_mesh_processing::detect_sharp_edges()`
 * @see `CGAL::Polygon_mesh_processing::connected_components()`
 * @see `CGAL::Polygon_mesh_processing::detect_vertex_incident_patches()`
 */
#ifdef DOXYGEN_RUNNING
template <typename PolygonMesh, typename FT,
          typename EdgeIsFeatureMap, typename PatchIdMap, typename NamedParameters>
#else
template <typename PolygonMesh,
          typename EdgeIsFeatureMap, typename PatchIdMap, typename NamedParameters = parameters::Default_named_parameters>
#endif
typename boost::graph_traits<PolygonMesh>::faces_size_type
sharp_edges_segmentation(const PolygonMesh& pmesh,
#ifdef DOXYGEN_RUNNING
                         FT angle_in_deg,
#else
                         typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT angle_in_deg,
#endif
                         EdgeIsFeatureMap edge_is_feature_map,
                         PatchIdMap patch_id_map,
                         const NamedParameters& np = parameters::default_values())
{
  detect_sharp_edges(pmesh, angle_in_deg, edge_is_feature_map, np);

  typename boost::graph_traits<PolygonMesh>::faces_size_type result =
    internal::detect_surface_patches(pmesh, patch_id_map, edge_is_feature_map, np);

  internal::vip_call(pmesh, patch_id_map,
    parameters::get_parameter(np, internal_np::vertex_incident_patches), edge_is_feature_map);

  return result;
}

} // end namespace PMP
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_H
