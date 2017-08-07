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
typename GT::Vector_3
facet_normal(const PolygonMesh& polygonMesh,
             const typename boost::graph_traits<PolygonMesh>::face_descriptor& f)
{
    return Polygon_mesh_processing::compute_face_normal(f,polygonMesh);
}


template <typename PolygonMesh, typename GT>
bool
is_sharp(PolygonMesh& polygonMesh,
         const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor& he,
         typename GT::FT cos_angle)
{
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    if(is_border(edge(he,polygonMesh),polygonMesh)){
        return false;
    }
    face_descriptor f1 = face(he,polygonMesh);
    face_descriptor f2 = face(opposite(he,polygonMesh),polygonMesh);

    const typename GT::Vector_3& n1 = facet_normal<PolygonMesh, GT>(polygonMesh, f1);
    const typename GT::Vector_3& n2 = facet_normal<PolygonMesh, GT>(polygonMesh, f2);

    if ( n1 * n2 <= cos_angle )
        return true;
    else
        return false;
}


template <typename PolygonMesh, typename PatchIdMap, typename EIF_map>
void
flood(PolygonMesh& polygonMesh,
      typename boost::graph_traits<PolygonMesh>::face_descriptor f,
      PatchIdMap& pid_map,
      const typename boost::property_traits<PatchIdMap>::value_type& patch_id,
      std::set<typename boost::graph_traits<PolygonMesh>::face_descriptor>& unsorted_faces,
      EIF_map& eif)
{
    typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor     face_descriptor;

    typedef std::set<face_descriptor> face_descriptor_set;
    typedef std::set<halfedge_descriptor> He_handle_set;
    // Initialize he_to_explore with halfedges of the starting facet
    He_handle_set he_to_explore;
    BOOST_FOREACH(halfedge_descriptor hd,
                  CGAL::halfedges_around_face(halfedge(f,polygonMesh), polygonMesh))
    {
        he_to_explore.insert(opposite(hd,polygonMesh));
    }

    // While there is something to explore
    while ( ! he_to_explore.empty() )
    {
        // Get next halfedge to explore
        halfedge_descriptor he = *(he_to_explore.begin());
        he_to_explore.erase(he_to_explore.begin());

        // If we don't go through a border of the patch
        if ( ! get(eif, edge(he, polygonMesh)) && ! is_border(he,polygonMesh))
        {
            face_descriptor explored_facet = face(he,polygonMesh);

            // Mark facet and delete it from unsorted
            put(pid_map, explored_facet, patch_id);
            unsorted_faces.erase(explored_facet);

            // Add/Remove facet's halfedge to/from explore list
            BOOST_FOREACH(halfedge_descriptor hd,
                          CGAL::halfedges_around_face(halfedge(explored_facet,polygonMesh),
                                                      polygonMesh))
            {
                halfedge_descriptor current_he = hd;

                // do not explore heh again
                if ( current_he == he ) { continue; }

                // if current_he is not in to_explore set, add it, otherwise remove it
                // (because we just explore the facet he_begin is pointing to)
                if ( he_to_explore.erase(current_he) == 0 )
                {
                    he_to_explore.insert(opposite(current_he,polygonMesh));
                }
            }
        }
    }
}

} //end internal

/*!
 * \ingroup PMP_detect_features_grp
 *
 * Detects the sharp edges of `p` by comparing with `angle_in_deg` and computes the number of sharp edges incident to each vertex.
 *
 * Property maps for CGAL::edge_is_feature_t and vertex_feature_degree_t should be either
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
 * \param angle_in_deg the floor angle.
 * \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 *    \cgalParamBegin{edge_is_feature_map}  a property map that will contain the sharp edges of `p` \cgalParamEnd
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
            = choose_param(get_param(np, internal_np::edge_is_feature),
                           get(CGAL::edge_is_feature, p));

    typedef typename boost::lookup_named_param_def <
            internal_np::vertex_feature_degree_t,
            NamedParameters,
            typename boost::property_map<PolygonMesh,vertex_feature_degree_t>::type//default
            > ::type                                               VNFE_map;
    VNFE_map vnfe
            = choose_param(get_param(np, internal_np::vertex_feature_degree),
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
 * based on the feature edges.
 *
 * A filled property map for CGAL::edge_is_feature_t should be either
 * available as an internal property map to `p` or provided as one of the Named Parameters.
 *
 * \tparam PolygonMesh a model of `FaceGraph`
 * \tparam PatchIdMap a model of `ReadWritePropertyMap` with
   `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type
   and the desired patch id, model of `CopyConstructible` as value type.
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param p the polygon mesh
 * \param patch_id_map a property map that will contain the surface patch ids for the faces of `p`.
 * \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 *    \cgalParamBegin{edge_is_feature_map}  a property map containing the sharp edges of `p` \cgalParamEnd
 *    \cgalParamBegin{maximum_number_of_patches} an `std::size_t`* the highest number of patches. This usually starts at 1. \cgalParamEnd
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
            = choose_param(get_param(np, internal_np::edge_is_feature),
                           get(CGAL::edge_is_feature, p));

    std::size_t dummy = 1;
    std::size_t* current_surface_index_ =
            boost::choose_param(get_param(np, internal_np::maximum_number_of_patches), &dummy);



    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
    typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
    typedef typename boost::property_traits<PatchIdMap>::value_type PatchId;
    // Initialize unsorted_faces
    typedef std::set<face_descriptor> face_descriptor_set;
    face_descriptor_set unsorted_faces;
    BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::face_descriptor fd, faces(p))
    {
        unsorted_faces.insert(fd);
    }

    // Flood
    while ( ! unsorted_faces.empty() )
    {
        face_descriptor f = *(unsorted_faces.begin());
        unsorted_faces.erase(unsorted_faces.begin());

        const PatchId patch_id = internal::generate_patch_id(PatchId(),
                                                             *current_surface_index_);
        put(patch_id_map, f, patch_id);
        internal::flood(p, f, patch_id_map, patch_id,unsorted_faces, eif);
        ++(*current_surface_index_);
    }
    return *current_surface_index_ - 1;
}


/*!
 * \ingroup PMP_detect_features_grp
 *
 * Gets for each vertex the patch ids of the incident faces.
 *
 * * A filled property map for CGAL::edge_is_feature_t should be either
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
 * \param patch_id_map a property map containing the surface patch ids for the faces of `p`. It must be already filled.
 * \param vertex_incident_patches_map a property map that will contain the patch ids of all the incident faces of each vertex of `p`.
 * \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 *    \cgalParamBegin{edge_is_feature_map}  a property map containing the sharp edges of `p` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * * @see `CGAL::Polygon_mesh_processing::detect_features()`
 */

template <typename PolygonMesh, typename PatchIdMap, typename VertexIncidentPatchesMap, typename NamedParameters>
void detect_incident_patches_of_vertices(PolygonMesh& p,
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
            = choose_param(get_param(np, internal_np::edge_is_feature),
                           get(CGAL::edge_is_feature, p));


    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

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
 * Detects the sharpd edges of `p` according to `angle_in_deg` and computes the corresponding
 * surface patch ids for each face.
 *
 * This function calls successively `CGAL::Polygon_mesh_processing::detect_sharp_edges()`,
 * `CGAL::Polygon_mesh_processing::detect_surface_patches()` and
 * `CGAL::Polygon_mesh_processing::detect_incident_patches_of_vertices()`
 *
 * Property maps for CGAL::edge_is_feature_t and vertex_feature_degree_t should be either
 * available as internal property maps to `p` or provided as Named Parameters.
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
 * \param angle_in_deg the floor angle.
 * \param patch_id_map a property map that will contain the surface patch ids for the faces of `p`.
 * \param vertex_incident_patches_map a property map that will contain the patch ids of all the incident faces of each vertex of `p`.
 * \param np optional \ref namedparameters described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
 *    \cgalParamBegin{edge_is_feature_map} a property map that will contain the sharp edges of `p` \cgalParamEnd
 *    \cgalParamBegin{vertex_feature_degree_map} a property map that will contain the number of adjacent feature edges for each vertex of `p` \cgalParamEnd
 *    \cgalParamBegin{maximum_number_of_patches} an `std::size_t` that that will contain the number of surface patches of `p`. \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @see `CGAL::Polygon_mesh_processing::detect_sharp_edges()`
 * @see `CGAL::Polygon_mesh_processing::detect_surface_patches()`
 * @see `CGAL::Polygon_mesh_processing::detect_incident_patches_of_vertices()`
 */
#ifdef DOXYGEN_RUNNING
template <typename PolygonMesh, typename FT, typename PatchIdMap, typename VertexIncidentPatchesMap, typename NamedParameters>
#else
template <typename PolygonMesh, typename PatchIdMap, typename VertexIncidentPatchesMap, typename NamedParameters>
#endif
void detect_features(PolygonMesh& p,
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
    detect_surface_patches(p, patch_id_map, np);
    detect_incident_patches_of_vertices(p, patch_id_map, vertex_incident_patches_map, np);
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
void detect_incident_patches_of_vertices(PolygonMesh& p,
                                      PatchIdMap& patch_id_map,
                                      VertexIncidentPatchesMap& vertex_incident_patches_map)
{
    detect_incident_patches_of_vertices(p, patch_id_map, vertex_incident_patches_map, parameters::all_default());
}
template <typename PolygonMesh, typename FT, typename PatchIdMap, typename VertexIncidentPatchesMap>
void detect_features(PolygonMesh& p,
                     FT angle_in_deg,
                     PatchIdMap& patch_id_map,
                     VertexIncidentPatchesMap& vertex_incident_patches_map)
{
    detect_features(p,angle_in_deg, patch_id_map, vertex_incident_patches_map, parameters::all_default());
}




} // end namespace PMP
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_DETECT_FEATURES_IN_POLYGON_MESH_H
