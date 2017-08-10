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
#include <CGAL/Mesh_3/properties.h>
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
template<typename PatchIdMap>
struct PatchIdMapWrapper
{
    typedef typename boost::property_traits<PatchIdMap>::category category;
    typedef typename boost::property_traits<PatchIdMap>::value_type value_type;
    typedef typename boost::property_traits<PatchIdMap>::reference reference;
    typedef typename boost::property_traits<PatchIdMap>::key_type key_type;

 PatchIdMap map;
 value_type offset;
 PatchIdMapWrapper(PatchIdMap& map, value_type offset)
     :map(map), offset(offset){}
};
template <typename PatchIdMap, typename Handle_type>
typename PatchIdMapWrapper<PatchIdMap>::value_type get(PatchIdMapWrapper<PatchIdMap>& map, Handle_type h)
{
    return get(map.map, h) - map.offset;
}

template <typename PatchIdMap, typename Handle_type>
void put(PatchIdMapWrapper<PatchIdMap>& map, Handle_type h,
         typename PatchIdMapWrapper<PatchIdMap>::value_type pid)
{
  put(map.map, h, pid + map.offset);
}

} //end internal

/*!
 * \ingroup PMP_detect_features_grp
 *
 * Detects the sharp edges of `p` according to `angle_in_deg` and computes the number of sharp edges incident to each vertex.
 *
 * Property maps for `CGAL::edge_is_feature_t` and `CGAL::vertex_feature_degree_t` should be either
 * available as internal property maps to `p` or provided as Named Parameters.
 *
 * \tparam PolygonMesh a model of `FaceGraph`
 * \tparam FT a number type. It is
 * either deduced from the `geom_traits` \ref namedparameters if provided,
 * or from the geometric traits class deduced from the point property map
 * of `PolygonMesh`.
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param p the polygon mesh
 * \param angle_in_deg the floor dihedral angle.
 * \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 *    \cgalParamBegin{edge_is_feature_map}  a property map that will contain the constrained-or-not status of each edge of `p` \cgalParamEnd
 *    \cgalParamBegin{vertex_feature_degree_map}  a property map that will contain the number of adjacent feature edges for each vertex of `p` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
#ifdef DOXYGEN_RUNNING
template <typename PolygonMesh, typename FT, typename NamedParameters>
#else
template <typename PolygonMesh, typename NamedParameters>
#endif
void detect_sharp_edges(PolygonMesh& p,
                        #ifdef DOXYGEN_RUNNING
                        FT angle_in_deg,
                        #else
                        typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT angle_in_deg,
                        #endif
                        const NamedParameters& np)
{
    //extract types from NPs
    typedef typename boost::lookup_named_param_def <
            internal_np::edge_is_feature_t,
            NamedParameters,
            typename boost::property_map<PolygonMesh,edge_is_feature_t>::type//default
            > ::type                                               EIF_map;
    EIF_map eif
            = boost::choose_param(get_param(np, internal_np::edge_is_feature),
                           get(CGAL::edge_is_feature, p));

    typedef typename boost::lookup_named_param_def <
            internal_np::vertex_feature_degree_t,
            NamedParameters,
            typename boost::property_map<PolygonMesh,vertex_feature_degree_t>::type//default
            > ::type                                               VNFE_map;
    VNFE_map vnfe
            = boost::choose_param(get_param(np, internal_np::vertex_feature_degree),
                           get(CGAL::vertex_feature_degree_t(), p));

    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
    typedef typename GT::FT FT;

    // Initialize vertices

    BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd, vertices(p))
    {
        put(vnfe,vd, 0);
    }

    FT cos_angle ( std::cos(CGAL::to_double(angle_in_deg) * CGAL_PI / 180.) );

    // Detect sharp edges
    BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::edge_descriptor ed, edges(p))
    {
        typename boost::graph_traits<PolygonMesh>::halfedge_descriptor he = halfedge(ed,p);
        if(is_border(he,p) || angle_in_deg == FT() ||
                (angle_in_deg != FT(180) && internal::is_sharp<PolygonMesh, GT>(p,he,cos_angle))
                )
        {
            put(eif, edge(he, p), true);

            put(vnfe, target(he,p), get(vnfe, target(he,p))+1);
            put(vnfe, source(he,p), get(vnfe, source(he,p))+1);
        }
    }
}

/*!
 * \ingroup PMP_detect_features_grp
 *
 * Computes for each face the index of the corresponding surface patch,
 * based on the feature edges which are considered as barriers between surface patches.
 * Filled property maps for `CGAL::edge_is_feature_t` and `CGAL::face_index_t`should be either
 * available as a internal property maps to `p` or provided as Named Parameters.
 *
 * \tparam PolygonMesh a model of `FaceGraph`
 * \tparam PatchIdMap a model of `ReadWritePropertyMap` with
   `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type
   and the desired patch id, model of `CopyConstructible` as value type.
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param p the polygon mesh
 * \param patch_id_map the property map that will contain the surface patch ids for the faces of `p`.
 * \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 *    \cgalParamBegin{edge_is_feature_map}  a property map containing the sharp edges of `p` \cgalParamEnd
 *    \cgalParamBegin{first_index} an std::size_t containing the index of the first surface patch of `p` \cgalParamEnd
 *    \cgalParamBegin{face_index_map} a property map containing the index of each face of `p` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \returns the number of surface patches.
 *
 * @see `CGAL::Polygon_mesh_processing::detect_features()`
 */
template <typename PolygonMesh, typename PatchIdMap, typename NamedParameters>
std::size_t
detect_surface_patches(PolygonMesh& p,
                            PatchIdMap& patch_id_map,
                            const NamedParameters& np)
{
    //extract types from NPs
    typedef typename boost::lookup_named_param_def <
            internal_np::edge_is_feature_t,
            NamedParameters,
            typename boost::property_map<PolygonMesh,edge_is_feature_t>::type//default
            > ::type                                               EIF_map;
    EIF_map eif
            = boost::choose_param(get_param(np, internal_np::edge_is_feature),
                           get(CGAL::edge_is_feature, p));
    typename GetFaceIndexMap<PolygonMesh, NamedParameters>::const_type
      fimap = boost::choose_param(get_param(np, internal_np::face_index),
                           get_const_property_map(boost::face_index, p));

    std::size_t offset = boost::choose_param(get_param(np, internal_np::first_index),
                                                      1);

    internal::PatchIdMapWrapper<PatchIdMap> wrapmap(patch_id_map, offset);
    return connected_components(p, wrapmap, parameters::edge_is_constrained_map(eif)
                                .face_index_map(fimap));

}


/*!
 * \ingroup PMP_detect_features_grp
 *
 * Collects the surface patches of the faces incident to each vertex
 *
 * * A filled property map for `CGAL::edge_is_feature_t` should be either
 * available as an internal property map to `p` or provided as one of the Named Parameters.
 *
 * \tparam PolygonMesh a model of `FaceGraph`
 * \tparam PatchIdMap a model of `ReadWritePropertyMap` with
   `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type
   and the desired patch id, model of `CopyConstructible` as value type.
 * \tparam VertexIncidentPatchesMap a model of `ReadWritePropertyMap` with
   `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type
   and a set of the desired patch id, model of `CopyConstructible` as value type.

 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param p the polygon mesh
 * \param patch_id_map the property map containing the surface patch ids for the faces of `p`. It must be already filled.
 * \param vertex_incident_patches_map a property map that will contain the patch ids of all the faces incident to each vertex of `p`.
 * \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 *    \cgalParamBegin{edge_is_feat  ure_map}  a property map containing the sharp edges of `p` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * * @see `CGAL::Polygon_mesh_processing::detect_features()`
 */

template <typename PolygonMesh, typename PatchIdMap, typename VertexIncidentPatchesMap, typename NamedParameters>
void detect_incident_patches(PolygonMesh& p,
                                      PatchIdMap& patch_id_map,
                                      VertexIncidentPatchesMap& vertex_incident_patches_map,
                                      const NamedParameters& np)
{
    //extract types from NPs
    typedef typename boost::lookup_named_param_def <
            internal_np::edge_is_feature_t,
            NamedParameters,
            typename boost::property_map<PolygonMesh,edge_is_feature_t>::type//default
            > ::type                                               EIF_map;
    EIF_map eif
            = boost::choose_param(get_param(np, internal_np::edge_is_feature),
                           get(CGAL::edge_is_feature, p));

    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor  vertex_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

    typedef typename boost::property_traits<PatchIdMap>::value_type PatchId;
    BOOST_FOREACH(vertex_descriptor vit,vertices(p))
    {
        // Look only at feature vertices
        if( ! get(eif, edge(halfedge(vit, p), p) )){ continue; }

        // Loop on incident facets of vit
        std::set<PatchId> set;
        BOOST_FOREACH(halfedge_descriptor he, halfedges_around_target(vit,p))
        {
            if( ! is_border(he,p) )
            {
                set.insert(get(patch_id_map,face(he,p)));
            }
            else if( ! is_border(opposite(he,p),p) )
            {
                set.insert(get(patch_id_map, face(opposite(he,p),p)));
            }
        }
        put(vertex_incident_patches_map, vit, set);
    }
}

/*!
 * \ingroup PMP_detect_features_grp
 *
 * Detects the sharp edges of `p` according to `angle_in_deg` and computes the corresponding
 * surface patch ids for each face.
 *
 * This function calls successively `CGAL::Polygon_mesh_processing::detect_sharp_edges()`,
 * `CGAL::Polygon_mesh_processing::detect_surface_patches()` and
 * `CGAL::Polygon_mesh_processing::detect_incident_patches()`
 *
 * Property maps for `CGAL::edge_is_feature_t`, `CGAL::vertex_feature_degree_t` and `CGAL::face_index_t`
 * should be either available as internal property maps to `p` or provided as Named Parameters.
 *
 * \tparam PolygonMesh a model of `FaceGraph`
 * \tparam FT a number type. It is
 * either deduced from the `geom_traits` \ref namedparameters if provided,
 * or from the geometric traits class deduced from the point property map
 * of `PolygonMesh`.
 * \tparam PatchIdMap a model of `ReadWritePropertyMap` with
   `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type
   and the desired patch id, model of `CopyConstructible` as value type.
 * \tparam VertexIncidentPatchesMap a model of `ReadWritePropertyMap` with
   `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type
   and a set of the desired patch id, model of `CopyConstructible` as value type.

 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param p the polygon mesh
 * \param angle_in_deg the floor dihedral angle.
 * \param patch_id_map the property map that will contain the surface patch ids for the faces of `p`.
 * \param vertex_incident_patches_map a property map that will contain the patch ids of all the faces incident to each vertex of `p`.
 * \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 *    \cgalParamBegin{edge_is_feature_map} a property map that will contain the constrained-or-not status of each edge of `p` \cgalParamEnd
 *    \cgalParamBegin{vertex_feature_degree_map} a property map that will contain the number of adjacent feature edges for each vertex of `p` \cgalParamEnd
 *    \cgalParamBegin{first_index} an std::size_t containing the index of the first surface patch of `p` \cgalParamEnd
 *    \cgalParamBegin{face_index_map} a property map containing the index of each face of `p` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \returns the number of surface patches.
 *
 * @see `CGAL::Polygon_mesh_processing::detect_sharp_edges()`
 * @see `CGAL::Polygon_mesh_processing::detect_surface_patches()`
 * @see `CGAL::Polygon_mesh_processing::detect_incident_patches()`
 */
#ifdef DOXYGEN_RUNNING
template <typename PolygonMesh, typename FT, typename PatchIdMap, typename VertexIncidentPatchesMap, typename NamedParameters>
#else
template <typename PolygonMesh, typename PatchIdMap, typename VertexIncidentPatchesMap, typename NamedParameters>
#endif
std::size_t detect_features(PolygonMesh& p,
                            #ifdef DOXYGEN_RUNNING
                            FT angle_in_deg,
                            #else
                            typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT angle_in_deg,
                            #endif
                            PatchIdMap& patch_id_map,
                            VertexIncidentPatchesMap& vertex_incident_patches_map,
                            const NamedParameters& np)
{
    detect_sharp_edges(p,angle_in_deg, np);
    std::size_t result = detect_surface_patches(p, patch_id_map, np);
    detect_incident_patches(p, patch_id_map, vertex_incident_patches_map, np);
    return result;
}

//Convenient overrides
template <typename PolygonMesh, typename FT>
void detect_sharp_edges(PolygonMesh& p,
                        FT angle_in_deg)
{
    detect_sharp_edges(p, angle_in_deg, parameters::all_default());
}

template <typename PolygonMesh, typename PatchIdMap>
typename boost::graph_traits<PolygonMesh>::faces_size_type
detect_surface_patches(PolygonMesh& p,
                            PatchIdMap& patch_id_map)
{
    return detect_surface_patches(p, patch_id_map, parameters::all_default());
}
template <typename PolygonMesh, typename PatchIdMap, typename VertexIncidentPatchesMap>
void detect_incident_patches(PolygonMesh& p,
                                      PatchIdMap& patch_id_map,
                                      VertexIncidentPatchesMap& vertex_incident_patches_map)
{
    detect_incident_patches(p, patch_id_map, vertex_incident_patches_map, parameters::all_default());
}
template <typename PolygonMesh, typename FT, typename PatchIdMap, typename VertexIncidentPatchesMap>
std::size_t detect_features(PolygonMesh& p,
                     FT angle_in_deg,
                     PatchIdMap& patch_id_map,
                     VertexIncidentPatchesMap& vertex_incident_patches_map)
{
    return detect_features(p,angle_in_deg, patch_id_map, vertex_incident_patches_map, parameters::all_default());
}




} // end namespace PMP
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_H
