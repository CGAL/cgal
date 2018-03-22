// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Laurent Rineau, Stephane Tayeb, Maxime Gimeno
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_H
#define CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_H

#include <CGAL/license/Polygon_mesh_processing/detect_features.h>

#include <CGAL/disable_warnings.h>

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

template <typename PolygonMesh, typename GT>
bool
is_sharp(PolygonMesh& polygonMesh,
         const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor& he,
         const typename GT::FT& cos_angle)
{
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  if(is_border(edge(he,polygonMesh),polygonMesh)){
    return false;
  }
  face_descriptor f1 = face(he,polygonMesh);
  face_descriptor f2 = face(opposite(he,polygonMesh),polygonMesh);

  const typename GT::Vector_3& n1 = Polygon_mesh_processing::compute_face_normal(f1,polygonMesh);
  const typename GT::Vector_3& n2 = Polygon_mesh_processing::compute_face_normal(f2,polygonMesh);

  if ( n1 * n2 <= cos_angle )
    return true;
  else
    return false;
}


//wrapper for patchid map.
template<typename PatchIdMap,
         typename ValueType = typename boost::property_traits<PatchIdMap>::value_type>
struct PatchIdMapWrapper
{
  typedef typename boost::property_traits<PatchIdMap>::category category;
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
  return get(map.map, h) - map.offset;
}

template <typename PatchIdMap, typename Handle_type, typename Int>
void put(PatchIdMapWrapper<PatchIdMap, Int>& map, Handle_type h,
         typename PatchIdMapWrapper<PatchIdMap, Int>::value_type pid)
{
  put(map.map, h, pid + map.offset);
}


//specialization for std::pair
template<typename PatchIdMap, typename Int>
struct PatchIdMapWrapper<PatchIdMap, std::pair<Int, Int> >
{
  typedef typename boost::property_traits<PatchIdMap>::category category;
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
  return get(map.map, h).first - map.offset;
}

template <typename PatchIdMap, typename Handle_type, typename Int>
void put(PatchIdMapWrapper<PatchIdMap, std::pair<Int, Int> >& map, Handle_type h,
         typename PatchIdMapWrapper<PatchIdMap, std::pair<Int, Int> >::value_type pid)
{
  put(map.map, h, std::pair<Int, Int>(pid+map.offset, 0));
}

template <typename PolygonMesh, typename PatchIdMap,
          typename EdgeIsFeatureMap, typename NamedParameters>
typename boost::graph_traits<PolygonMesh>::faces_size_type
detect_surface_patches(PolygonMesh& p,
                       PatchIdMap patch_id_map,
                       EdgeIsFeatureMap eif,
                       const NamedParameters& np)
{
  //extract types from NPs
  typename GetFaceIndexMap<PolygonMesh, NamedParameters>::const_type
          fimap = boost::choose_param(get_param(np, internal_np::face_index),
                                      get_const_property_map(boost::face_index, p));

  int offset = static_cast<int>(
          boost::choose_param(get_param(np, internal_np::first_index),
          1));

  internal::PatchIdMapWrapper<PatchIdMap,
                              typename boost::property_traits<PatchIdMap>::value_type>
          wrapmap(patch_id_map, offset);
  return connected_components(p, wrapmap,
                              parameters::edge_is_constrained_map(eif)
                             .face_index_map(fimap));

}
template <typename PolygonMesh, typename EdgeIsFeatureMap, typename PatchIdMap>
typename boost::graph_traits<PolygonMesh>::faces_size_type
detect_surface_patches(PolygonMesh& p,
                       PatchIdMap patch_id_map,
                       EdgeIsFeatureMap eif)
{
  return detect_surface_patches(p, patch_id_map, eif, parameters::all_default());
}


template<typename GT,
         typename FT,
         typename PolygonMesh,
         typename EIFMap,
         typename VNFEMap>
 void sharp_call(PolygonMesh& pmesh,
                 FT angle_in_deg,
                 EIFMap edge_is_feature_map,
                 VNFEMap vnfe)
{
  // Initialize vertices
  BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd,
                vertices(pmesh))
  {
    put(vnfe, vd, 0);
  }
  FT cos_angle ( std::cos(CGAL::to_double(angle_in_deg) * CGAL_PI / 180.) );

  // Detect sharp edges
  BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::edge_descriptor ed, edges(pmesh))
  {
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor he = halfedge(ed,pmesh);
    if(is_border_edge(he,pmesh)
      || angle_in_deg == FT()
      || (angle_in_deg != FT(180) && internal::is_sharp<PolygonMesh, GT>(pmesh,he,cos_angle))
      )
    {
      put(edge_is_feature_map, edge(he, pmesh), true);
      put(vnfe, target(he,pmesh), get(vnfe, target(he,pmesh))+1);
      put(vnfe, source(he,pmesh), get(vnfe, source(he,pmesh))+1);
    }
  }
}


template<typename GT,
         typename FT,
         typename PolygonMesh,
         typename EIFMap>
 void sharp_call(PolygonMesh& pmesh,
                 FT& angle_in_deg,
                 EIFMap edge_is_feature_map,
                 const boost::param_not_found&)
{
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor     edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  FT cos_angle ( std::cos(CGAL::to_double(angle_in_deg) * CGAL_PI / 180.) );

  // Detect sharp edges
  BOOST_FOREACH(edge_descriptor ed, edges(pmesh))
  {
    halfedge_descriptor he = halfedge(ed,pmesh);
    if(is_border_edge(he,pmesh)
      || angle_in_deg == FT()
      || (angle_in_deg != FT(180) && internal::is_sharp<PolygonMesh, GT>(pmesh,he,cos_angle))
      )
    {
      put(edge_is_feature_map, edge(he, pmesh), true);
    }
  }
}
} //end internal

/*!
 * \ingroup PMP_detect_features_grp
 *
 * detects and marks the edges that are considered to be sharp with respect to the given angle bound.
 * `angle_in_deg` gives the maximum angle (in degrees) between the two normal vectors of adjacent triangles.
 * For an edge of the input polygon mesh, if the angle between the two normal vectors of its incident facets is bigger
 * than the given bound, then the edge is marked as being a feature edge.
 *
 * Also computes the number of sharp edges incident to each vertex, if `vertex_feature_degree_map` is provided.
 *
 * \tparam PolygonMesh a model of `HalfedgeListGraph`
 * \tparam FT a number type. It is
 * either deduced from the `geom_traits` \ref pmp_namedparameters "Named Parameters" if provided,
 * or from the geometric traits class deduced from the point property map
 * of `PolygonMesh`.
 * \tparam EdgeIsFeatureMap a model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
 *  as key type and `bool` as value type. It should be default constructible.
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param pmesh the polygon mesh
 * \param angle_in_deg the dihedral angle bound
 * \param edge_is_feature_map the property map that will contain the sharp-or-not status of each edge of `pmesh`
 * \param np optional \ref pmp_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 *    \cgalParamBegin{vertex_feature_degree_map} a property map that will contain the number of adjacent feature edges for
 *  each vertex of `pmesh` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
#ifdef DOXYGEN_RUNNING
template <typename PolygonMesh, typename FT,
          typename EdgeIsFeatureMap, typename NamedParameters>
#else
template <typename PolygonMesh, typename EdgeIsFeatureMap, typename NamedParameters>
#endif
void detect_sharp_edges(PolygonMesh& pmesh,
#ifdef DOXYGEN_RUNNING
    FT angle_in_deg,
#else
    typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT angle_in_deg,
#endif
    EdgeIsFeatureMap edge_is_feature_map,
    const NamedParameters& np)
{
  //extract types from NPs
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
  typedef typename GetGeomTraits<PolygonMesh, GT>::type::FT          FT;

  internal::sharp_call<GT, FT>(pmesh, angle_in_deg, edge_is_feature_map,
                               get_param(np, internal_np::vertex_feature_degree));
}


/*!
 * \ingroup PMP_detect_features_grp
 *
 * collects the surface patches of the faces incident to each vertex of the input polygon mesh.
 *
 * \tparam PolygonMesh a model of `HalfedgeListGraph`
 * \tparam PatchIdMap a model of `ReadPropertyMap` with
   `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type
   and the desired patch id, model of `CopyConstructible` as value type.
 * \tparam VertexIncidentPatchesMap a model of mutable `LvaluePropertyMap` with
   `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type. Its value type
   must be a container of `boost::property_traits<PatchIdMap>::%value_type` and have a function `insert()`.
   A `std::set` or a `boost::unordered_set` are recommended, as a patch index may be
   inserted several times.
 * \tparam EdgeIsFeatureMap a model of `ReadPropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
 *  as key type and `bool` as value type.
 * \param pmesh the polygon mesh
 * \param patch_id_map the property map containing the surface patch ids for the faces of `pmesh`. It must be already filled.
 * \param vertex_incident_patches_map a property map that will contain the patch ids of all the faces incident to each vertex of `pmesh`.
 * \param edge_is_feature_map a filled property map that will contain the sharp-or-not status of each edge of `pmesh`
 *
 * @see `CGAL::Polygon_mesh_processing::sharp_edges_segmentation()`
 */

template <typename PolygonMesh, typename PatchIdMap,
          typename VertexIncidentPatchesMap, typename EdgeIsFeatureMap>
void detect_vertex_incident_patches(PolygonMesh& pmesh,
                             const PatchIdMap patch_id_map,
                             VertexIncidentPatchesMap vertex_incident_patches_map,
                             const EdgeIsFeatureMap edge_is_feature_map)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  BOOST_FOREACH(vertex_descriptor vit,vertices(pmesh))
  {
    // Look only at feature vertices
    if( ! get(edge_is_feature_map, edge(halfedge(vit, pmesh), pmesh) ))
    { continue; }

    // Loop on incident facets of vit
    typename VertexIncidentPatchesMap::value_type&
      id_set = vertex_incident_patches_map[vit];
    BOOST_FOREACH(halfedge_descriptor he, halfedges_around_target(vit,pmesh))
    {
      if( ! is_border(he,pmesh) )
      {
        id_set.insert(get(patch_id_map, face(he, pmesh)));
      }
      else if( ! is_border(opposite(he,pmesh),pmesh) )
      {
        id_set.insert(get(patch_id_map, face(opposite(he, pmesh), pmesh)));
      }
    }
  }
}


namespace internal
{
  template<typename PolygonMesh,
              typename PIDMap,
              typename VIPMap,
              typename EIFMap>
  void vip_call(PolygonMesh& mesh, PIDMap pid, VIPMap vip, EIFMap eif )
  {
    CGAL::Polygon_mesh_processing::detect_vertex_incident_patches(mesh, pid, vip, eif);
  }

  template<typename PolygonMesh,
           typename PIDMap,
           typename EIFMap>
  void vip_call(PolygonMesh&, PIDMap, const boost::param_not_found&, EIFMap)
  {
    //do nothing when the parameter is not given
  }
}//end internal
/*!
 * \ingroup PMP_detect_features_grp
 *
 * This function calls successively `CGAL::Polygon_mesh_processing::detect_sharp_edges()`,
 * `CGAL::Polygon_mesh_processing::connected_components()` and
 * `CGAL::Polygon_mesh_processing::detect_vertex_incident_patches()`
 *
 * It detects and marks the sharp edges of `pmesh` according to `angle_in_deg`.
 * The sharp edges define a segmentation of `pmesh`, that is done by
 * computing a
 * surface patch id for each face.
 *
 * A property map for `CGAL::face_index_t`should be either available
 * as an internal property map to `pmesh` or provided as one of the Named Parameters.
 *
 * \tparam PolygonMesh a model of `FaceGraph`
 * \tparam FT a number type. It is
 * either deduced from the `geom_traits` \ref pmp_namedparameters "Named Parameters" if provided,
 * or from the geometric traits class deduced from the point property map
 * of `PolygonMesh`.
 * \tparam EdgeIsFeatureMap a model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
 * \tparam PatchIdMap a model of `ReadWritePropertyMap` with
   `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type
   and the desired patch id, model of `CopyConstructible` as value type.
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param pmesh the polygon mesh
 * \param angle_in_deg the dihedral angle bound
 * \param edge_is_feature_map the property map that will contain the sharp-or-not status of each edge of `pmesh`
 * \param patch_id_map the property map that will contain the surface patch ids for the faces of `pmesh`.
 * \param np optional \ref pmp_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 *    \cgalParamBegin{vertex_feature_degree_map} a property map that will contain the number of adjacent feature edges for each vertex of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{first_index} a `std::size_t` containing the index of the first surface patch of `pmesh`.
 *      The patches will be numbered on [first_index; first_index + num_patches], where num_patches is the number of surface patches \cgalParamEnd
 *    \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{vertex_incident_patches_map} a property map that will contain the patch ids of all the faces incident to each vertex of `pmesh`. \cgalParamEnd
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
          typename EdgeIsFeatureMap, typename PatchIdMap, typename NamedParameters>
#endif
typename boost::graph_traits<PolygonMesh>::faces_size_type
sharp_edges_segmentation(PolygonMesh& pmesh,
#ifdef DOXYGEN_RUNNING
      FT angle_in_deg,
#else
      typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT angle_in_deg,
#endif
      EdgeIsFeatureMap edge_is_feature_map,
      PatchIdMap patch_id_map,
      const NamedParameters& np)
{
    detect_sharp_edges(pmesh, angle_in_deg, edge_is_feature_map, np);

    typename boost::graph_traits<PolygonMesh>::faces_size_type result =
      internal::detect_surface_patches(pmesh, patch_id_map, edge_is_feature_map, np);

    internal::vip_call(pmesh, patch_id_map,
      get_param(np, internal_np::vertex_incident_patches), edge_is_feature_map);

    return result;
}

//Convenient overrides
template <typename PolygonMesh, typename EdgeIsFeatureMap, typename FT>
void detect_sharp_edges(PolygonMesh& p,
                        FT angle_in_deg,
                        EdgeIsFeatureMap edge_is_feature_map)
{
  detect_sharp_edges(p, angle_in_deg, edge_is_feature_map,
                     parameters::all_default());
}

template <typename PolygonMesh, typename FT,
          typename EdgeIsFeatureMap, typename PatchIdMap>
typename boost::graph_traits<PolygonMesh>::faces_size_type
sharp_edges_segmentation(PolygonMesh& p,
                         FT angle_in_deg,
                         EdgeIsFeatureMap edge_is_feature_map,
                         PatchIdMap patch_id_map)
{
  return sharp_edges_segmentation(p, angle_in_deg, edge_is_feature_map, patch_id_map,
                                 parameters::all_default());
}




} // end namespace PMP
} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_H
