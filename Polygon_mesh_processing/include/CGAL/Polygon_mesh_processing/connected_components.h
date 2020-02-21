// Copyright (c) 2011, 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot and Andreas Fabri

#ifndef CGAL_POLYGON_MESH_PROCESSING_CONNECTED_COMPONENTS_H
#define CGAL_POLYGON_MESH_PROCESSING_CONNECTED_COMPONENTS_H

#include <CGAL/license/Polygon_mesh_processing/connected_components.h>

#include <CGAL/disable_warnings.h>

#include<set>
#include<vector>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/vector_property_map.hpp>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/Dual.h>
#include <CGAL/Default.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/iterator.h>
#include <CGAL/tuple.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif


namespace CGAL {
namespace Polygon_mesh_processing{
namespace internal {

  struct MoreSecond
  {
    template <typename T1, typename T2>
    bool operator()(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) const {
      return a.second > b.second;
    }
  };

    // A property map 
    template <typename G>
    struct No_constraint {
      friend bool get(No_constraint<G>, typename boost::graph_traits<G>::edge_descriptor)
      {
        return false;
      }
    };

    // A functor
    template <typename G, typename EdgeConstraintMap = No_constraint<G> >
    struct No_border {
      No_border()
      {}

      No_border(const G & g, EdgeConstraintMap ecm = EdgeConstraintMap())
        : g(&g), ecm(ecm)
      {}

      bool operator()(typename boost::graph_traits<G>::edge_descriptor e) const {
        if (!is_border(e, *g)){
          return !get(ecm, e);
        }
        return false;
      }

      const G* g;
      EdgeConstraintMap ecm;
    };

} // namespace internal

/*!
 * \ingroup keep_connected_components_grp
 *  discovers all the faces in the same connected component as `seed_face` and records them in `out`.
 * `seed_face` will also be added in `out`.
 *
 *  \tparam PolygonMesh a model of `FaceGraph`
 *  \tparam FaceOutputIterator a model of `OutputIterator` that accepts
        faces of type
        `boost::graph_traits<PolygonMesh>::%face_descriptor`.
 *  \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 *  \param seed_face a face of `pmesh` from which exploration starts to detect the connected component
           that contains it
 *  \param pmesh the polygon mesh
 *  \param out the output iterator that collects faces from the same connected component as `seed_face`
 *  \param np optional \ref pmp_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{edge_is_constrained_map}  a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 *  \returns the output iterator.
 *
 */
template <typename PolygonMesh
          , typename FaceOutputIterator
          , typename NamedParameters
          >
FaceOutputIterator
connected_component(typename boost::graph_traits<PolygonMesh>::face_descriptor seed_face
                    , const PolygonMesh& pmesh
                    , FaceOutputIterator out
                    , const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters,
    internal::No_constraint<PolygonMesh>//default
  > ::type                                               EdgeConstraintMap;
  EdgeConstraintMap ecmap
    = choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
                   internal::No_constraint<PolygonMesh>());

  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  std::set<face_descriptor> already_processed;
  std::vector< face_descriptor > stack;
  stack.push_back(seed_face);
  while (!stack.empty())
    {
      seed_face=stack.back();
      stack.pop_back();
      if (!already_processed.insert(seed_face).second) continue;
      *out++=seed_face;
      for(halfedge_descriptor hd :
                    halfedges_around_face(halfedge(seed_face, pmesh), pmesh) )
      {
        if(! get(ecmap, edge(hd, pmesh))){
          face_descriptor neighbor = face( opposite(hd, pmesh), pmesh );
          if ( neighbor != boost::graph_traits<PolygonMesh>::null_face() )
            stack.push_back(neighbor);
        }
      }
    }
  return out;
}

template <typename PolygonMesh, typename OutputIterator>
OutputIterator
connected_component(typename boost::graph_traits<PolygonMesh>::face_descriptor seed_face,
                    const PolygonMesh& pmesh,
                    OutputIterator out)
{
  return connected_component(seed_face, pmesh, out,
          CGAL::Polygon_mesh_processing::parameters::all_default());
}

/*!
 * \ingroup keep_connected_components_grp
 *  computes for each face the index of the corresponding connected component.
 *
 * If `PolygonMesh`
 *  has an internal non writable property map
 *  for `CGAL::face_index_t` and no `face_index_map` is given
 *  as a named parameter, then the internal one must be initialized; otherwise, it will be.
 *
 *  \tparam PolygonMesh a model of `FaceListGraph`
 *  \tparam FaceComponentMap a model of `WritablePropertyMap` with
        `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
        `boost::graph_traits<PolygonMesh>::%faces_size_type` as value type.
 *  \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"

 * \param pmesh the polygon mesh
 * \param fcm the property map with indices of components associated to faces in `pmesh`
 * \param np optional \ref pmp_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{edge_is_constrained_map} a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 *  \returns the number of connected components.
 */

template <typename PolygonMesh
        , typename FaceComponentMap
        , typename NamedParameters
>
typename boost::property_traits<FaceComponentMap>::value_type
connected_components(const PolygonMesh& pmesh,
                     FaceComponentMap fcm,
                     const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters,
    internal::No_constraint<PolygonMesh>//default
  > ::type                                               EdgeConstraintMap;
  
  EdgeConstraintMap ecmap
    = choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
                   internal::No_constraint<PolygonMesh>());

  typedef typename Default_face_index_map<NamedParameters, PolygonMesh>::type FaceIndexMap;
  
  FaceIndexMap  fimap = get_initialized_face_index_map(pmesh, np);
  
  typename boost::property_traits<FaceComponentMap>::value_type i=0;
  std::vector<bool> handled(num_faces(pmesh), false);
  for (face_descriptor f : faces(pmesh))
  {
    if (handled[get(fimap,f)]) continue;
    std::vector<face_descriptor> queue;
    queue.push_back(f);
    while(!queue.empty())
    {
      face_descriptor fq = queue.back();
      queue.pop_back();
      typename boost::property_traits<FaceIndexMap>::value_type  fq_id = get(fimap,fq);
      if ( handled[fq_id]) continue;
      handled[fq_id]=true;
      put(fcm, fq, i);
      for (halfedge_descriptor h : halfedges_around_face(halfedge(fq, pmesh), pmesh))
      {
        if ( get(ecmap, edge(h, pmesh)) ) continue;
        halfedge_descriptor opp = opposite(h, pmesh);
        face_descriptor fqo = face(opp, pmesh);
        if ( fqo != GT::null_face() )
        {
          if ( !handled[get(fimap,fqo)] )
            queue.push_back(fqo);
        }
      }
    }
    ++i;
  }
  return i;
}

template <typename PolygonMesh, typename FaceComponentMap>
typename boost::property_traits<FaceComponentMap>::value_type
connected_components(const PolygonMesh& pmesh,
                     FaceComponentMap fcm)
{
  return CGAL::Polygon_mesh_processing::connected_components(pmesh, fcm, CGAL::parameters::all_default());
}

template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap
        , typename NamedParameters>
void keep_connected_components(PolygonMesh& pmesh
                              , const ComponentRange& components_to_keep
                              , const FaceComponentMap& fcm
                              , const NamedParameters& np);

namespace internal {

//  /*!
//  * \ingroup keep_connected_components_grp
//  *  returns the number of connected components in the mesh.
//  *
//  *  A property map for `CGAL::face_index_t` must be either available as an internal property map
//  *  to `pmesh` or provided as one of the \ref pmp_namedparameters "Named Parameters".
//  *
//  *  \tparam PolygonMesh a model of `FaceGraph`
//  *  \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
//  *
//  *  \param pmesh the polygon mesh
//  *  \param np optional \ref pmp_namedparameters "Named Parameters" described below
//  *
//  * \cgalNamedParamsBegin
//  *  \cgalParamBegin{edge_is_constrained_map}  a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
//  *  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
//  * \cgalNamedParamsEnd
//  *
//  * \returns the output iterator.
//  *
//  */
template <typename PolygonMesh,
          typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
std::size_t number_of_connected_components(const PolygonMesh& pmesh,
                                           const CGAL_PMP_NP_CLASS& np)
{
  typedef CGAL::dynamic_face_property_t<std::size_t>                                Face_property_tag;
  typedef typename boost::property_map<PolygonMesh, Face_property_tag >::const_type Patch_ids_map;

  Patch_ids_map patch_ids_map = get(Face_property_tag(), pmesh);

  return CGAL::Polygon_mesh_processing::connected_components(pmesh, patch_ids_map, np);
}

template <typename PolygonMesh>
std::size_t number_of_connected_components(const PolygonMesh& pmesh)
{
  return internal::number_of_connected_components(pmesh, CGAL::parameters::all_default());
}

} // end namespace internal

/*!
 * \ingroup keep_connected_components_grp
 *
 * removes the small connected components and all isolated vertices.
 * Keep the `nb_components_to_keep` largest connected components, where the size of a connected
 * component is computed as the sum of the individual sizes of all the faces of the connected component.
 * By default, the size of a face is `1` (and thus the size of a connected component is the number
 * of faces it contains), but it is also possible to pass custom sizes, such as the area of the face.
 *
* If `PolygonMesh`
*  has a non writable internal property map
*  for `CGAL::face_index_t` or `CGAL::vertex_index_t` and no `face_index_map` (respectively `vertex_index_map`) is given
*  as a named parameter, then the internal one(s) must be initialized. Otherwise, it will be.
 *
 * \tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param pmesh the polygon mesh
 * \param nb_components_to_keep the number of components to be kept. If this number is larger than
 *                              the number of components in the mesh, all components are kept.
 * \param np optional \ref pmp_namedparameters "Named Parameters", amongst those described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{edge_is_constrained_map} a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{face_size_map}
 *      a property map containing a size for each face of `pmesh`. The value type of this property map
 *      is chosen by the user, but must be constructible from `0` and support `operator+=()` and
 *      comparisons.
 *    \cgalParamEnd
 *    \cgalParamBegin{dry_run}
 *      a Boolean parameter. If set to `true`, the mesh will not be altered, but the number
 *      of components that would be removed is returned. The default value is `false`.
 *    \cgalParamEnd
 *    \cgalParamBegin{output_iterator} a model of `OutputIterator` with value type `face_descriptor`.
 *      When using the "dry run" mode (see parameter `dry_run`), faces that would be removed by the
 *      algorithm can be collected with this output iterator.
 *    \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 *  \return the number of connected components removed (ignoring isolated vertices).
 */
template <typename PolygonMesh,
          typename NamedParameters>
std::size_t keep_largest_connected_components(PolygonMesh& pmesh,
                                              std::size_t nb_components_to_keep,
                                              const NamedParameters& np)
{
  typedef PolygonMesh                                                   PM;
  typedef typename boost::graph_traits<PM>::face_descriptor             face_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename Default_face_index_map<NamedParameters, PolygonMesh>::type FaceIndexMap;

  FaceIndexMap fimap = get_initialized_face_index_map(pmesh, np);
  // FaceSizeMap
  typedef typename internal_np::Lookup_named_param_def<internal_np::face_size_map_t,
                                                 NamedParameters,
                                                 Constant_property_map<face_descriptor, std::size_t> // default
                                                >::type                  FaceSizeMap;
  typedef typename boost::property_traits<FaceSizeMap>::value_type       Face_size;

  FaceSizeMap face_size_pmap = choose_parameter(get_parameter(np, internal_np::face_size_map),
                                                Constant_property_map<face_descriptor, std::size_t>(1));

  const bool dry_run = choose_parameter(get_parameter(np, internal_np::dry_run), false);

  typedef typename internal_np::Lookup_named_param_def<internal_np::output_iterator_t,
                                                       NamedParameters,
                                                       Emptyset_iterator>::type Output_iterator;
  Output_iterator out = choose_parameter(get_parameter(np, internal_np::output_iterator), Emptyset_iterator());

  // vector_property_map
  boost::vector_property_map<std::size_t, FaceIndexMap> face_cc(fimap);
  std::size_t num = connected_components(pmesh, face_cc, np);

  // Even if we do not want to keep anything we need to first
  // calculate the number of existing connected_components to get the
  // correct return value.
  if(nb_components_to_keep == 0)
  {
    CGAL::clear(pmesh);
    return num;
  }

  if(nb_components_to_keep >= num)
    return 0;

  std::vector<std::pair<std::size_t, Face_size> > component_size(num);

  for(std::size_t i=0; i < num; i++)
    component_size[i] = std::make_pair(i, Face_size(0));

  for(face_descriptor f : faces(pmesh))
    component_size[face_cc[f]].second += get(face_size_pmap, f);

  // we sort the range [0, num) by component size
  std::sort(component_size.begin(), component_size.end(), internal::MoreSecond());

  if(dry_run)
  {
    std::vector<bool> is_to_be_removed(num, false);
    for(std::size_t i=0; i<nb_components_to_keep; ++i)
      is_to_be_removed[component_size[i].first] = true;

    for(face_descriptor f : faces(pmesh))
      if(is_to_be_removed[face_cc[f]])
        *out++ = f;
  }
  else
  {
    std::vector<std::size_t> cc_to_keep;
    for(std::size_t i=0; i<nb_components_to_keep; ++i)
      cc_to_keep.push_back(component_size[i].first);

    keep_connected_components(pmesh, cc_to_keep, face_cc, np);
  }

  return num - nb_components_to_keep;
}

template <typename PolygonMesh>
std::size_t keep_largest_connected_components(PolygonMesh& pmesh,
                                              std::size_t nb_components_to_keep)
{
  return keep_largest_connected_components(pmesh,
    nb_components_to_keep,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

/*!
 * \ingroup keep_connected_components_grp
 * removes connected components whose size is (strictly) smaller than a given threshold value,
 * where the size of a connected component is computed as the sum of the individual sizes
 * of all the faces of the connected component. By default, the size of a face is `1` (and thus
 * the size of a connected component is the number of faces it contains), but it is also possible
 * to pass custom sizes, such as the area of the face.
 *
* If `PolygonMesh`
*  has a non writable internal property map
*  for `CGAL::face_index_t` or `CGAL::vertex_index_t` and no `face_index_map` (respectively `vertex_index_map`) is given
*  as a named parameter, then the internal one(s) must be initialized. Otherwise, it will be.
 *
 * \tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
 * \tparam ThresholdValueType the type of the threshold value
 * \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * \param pmesh the polygon mesh
 * \param threshold_value any connected component with a size (strictly) smaller than this value will be discarded
 * \param np optional \ref pmp_namedparameters "Named Parameters", amongst those described below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{edge_is_constrained_map} a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{face_size_map}
 *      a property map containing a size for each face of `pmesh`. The value type of this property map
 *      is chosen by the user, but must be constructible from `0` and support `operator+=()` and
 *      comparisons.
 *    \cgalParamEnd
 *    \cgalParamBegin{dry_run}
 *      a Boolean parameter. If set to `true`, the mesh will not be altered, but the number
 *      of components that would be removed is returned. The default value is `false`.
 *    \cgalParamEnd
 *    \cgalParamBegin{output_iterator} a model of `OutputIterator` with value type `face_descriptor`.
 *      When using the "dry run" mode (see parameter `dry_run`), faces that would be removed by the
 *      algorithm can be collected with this output iterator.
 *    \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \pre If a face size property map is passed by the user, `ThresholdValueType` must be the same
 *      type as the value type of the property map. Otherwise, `ThresholdValueType` must be `std::size_t`.
 *
 *  \return the number of connected components removed (ignoring isolated vertices).
 */
template <typename PolygonMesh,
          typename ThresholdValueType,
          typename NamedParameters>
std::size_t keep_large_connected_components(PolygonMesh& pmesh,
                                            const ThresholdValueType threshold_value,
                                            const NamedParameters& np)
{
  typedef PolygonMesh                                                     PM;
  typedef typename boost::graph_traits<PM>::face_descriptor               face_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename Default_face_index_map<NamedParameters, PolygonMesh>::type FaceIndexMap;
  FaceIndexMap fim =
      CGAL::get_initialized_face_index_map(pmesh, np);

  typedef typename internal_np::Lookup_named_param_def<internal_np::face_size_map_t,
                                                 NamedParameters,
                                                 Constant_property_map<face_descriptor, std::size_t> // default
                                                >::type                   FaceSizeMap;
  typedef typename boost::property_traits<FaceSizeMap>::value_type        Face_size;

  CGAL_static_assertion((std::is_convertible<ThresholdValueType, Face_size>::value));

  typedef typename internal_np::Lookup_named_param_def<internal_np::output_iterator_t,
                                                       NamedParameters,
                                                       Emptyset_iterator>::type Output_iterator;

  FaceSizeMap face_size_pmap = choose_parameter(get_parameter(np, internal_np::face_size_map),
                                           Constant_property_map<face_descriptor, std::size_t>(1));
  const bool dry_run = choose_parameter(get_parameter(np, internal_np::dry_run), false);
  Output_iterator out = choose_parameter(get_parameter(np, internal_np::output_iterator), Emptyset_iterator());

  // vector_property_map
  boost::vector_property_map<std::size_t, FaceIndexMap> face_cc(fim);
  std::size_t num = connected_components(pmesh, face_cc, np);
  std::vector<Face_size> component_size(num, 0);

  for(face_descriptor f : faces(pmesh))
    component_size[face_cc[f]] += get(face_size_pmap, f);

  const Face_size thresh = threshold_value;
  std::vector<bool> is_to_be_kept(num, false);
  std::size_t res = 0;

  for(std::size_t i=0; i<num; ++i)
  {
    if(component_size[i] >= thresh)
    {
      is_to_be_kept[i] = true;
      ++res;
    }
  }

  if(dry_run)
  {
    for(face_descriptor f : faces(pmesh))
      if(!is_to_be_kept[face_cc[f]])
        *out++ = f;
  }
  else
  {
    std::vector<std::size_t> ccs_to_keep;
    for(std::size_t i=0; i<num; ++i)
      if(is_to_be_kept[i])
        ccs_to_keep.push_back(i);

    keep_connected_components(pmesh, ccs_to_keep, face_cc, np);
  }

  return num - res;
}

template <typename PolygonMesh>
std::size_t keep_large_connected_components(PolygonMesh& pmesh,
                                            std::size_t threshold_components_to_keep)
{
  return keep_large_connected_components(pmesh,
    threshold_components_to_keep,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}


template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap
        , typename NamedParameters>
void keep_or_remove_connected_components(PolygonMesh& pmesh
                                        , const ComponentRange& components_to_keep
                                        , const FaceComponentMap& fcm
                                        , bool  keep
                                        , const NamedParameters& np)
{
  typedef PolygonMesh PM;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor   face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_iterator     face_iterator;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_iterator   vertex_iterator;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor   edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_iterator     edge_iterator;

  //VertexIndexMap
  typedef typename Default_vertex_index_map<NamedParameters, PM>::type VertexIndexMap;
   VertexIndexMap vim = get_initialized_vertex_index_map(pmesh, np);

  std::set<std::size_t> cc_to_keep;
  for(std::size_t i : components_to_keep)
    cc_to_keep.insert(i);

  boost::vector_property_map<bool, VertexIndexMap> keep_vertex(vim);
  for(vertex_descriptor v : vertices(pmesh)){
    keep_vertex[v] = false;
  }
  for(face_descriptor f : faces(pmesh)){
    if (cc_to_keep.find(get(fcm,f)) != cc_to_keep.end())
      put(fcm, f, keep ? 1 : 0);
    else
      put(fcm, f, keep ? 0 : 1);
  }

  for(face_descriptor f : faces(pmesh)){
    if (get(fcm, f) == 1){
      for(halfedge_descriptor h : halfedges_around_face(halfedge(f, pmesh), pmesh)){
        vertex_descriptor v = target(h, pmesh);
        keep_vertex[v] = true;
      }
    }
  }

  edge_iterator eb, ee;
  for (boost::tie(eb, ee) = edges(pmesh); eb != ee;)
  {
    edge_descriptor e = *eb;
    ++eb;
    vertex_descriptor v = source(e, pmesh);
    vertex_descriptor w = target(e, pmesh);
    halfedge_descriptor h = halfedge(e, pmesh);
    halfedge_descriptor oh = opposite(h, pmesh);
    if (!keep_vertex[v] && !keep_vertex[w]){
      // don't care about connectivity
      // As vertices are not kept the faces and vertices will be removed later
      remove_edge(e, pmesh);
    }
    else if (keep_vertex[v] && keep_vertex[w]){
      face_descriptor fh = face(h, pmesh), ofh = face(oh, pmesh);
      if (is_border(h, pmesh) && is_border(oh, pmesh)){
#ifdef CGAL_CC_DEBUG
        std::cerr << "null_face on both sides of " << e << " is kept\n";
#endif
      }
      else if ((is_border(oh, pmesh) && get(fcm,fh)) ||
        (is_border(h, pmesh) && get(fcm,ofh)) ||
        (!is_border(oh, pmesh) && !is_border(h, pmesh) && get(fcm,fh) && get(fcm,ofh))){
        // do nothing
      }
      else if (!is_border(h, pmesh) && get(fcm,fh) && !is_border(oh, pmesh) && !get(fcm,ofh)){
        set_face(oh, boost::graph_traits<PolygonMesh>::null_face(), pmesh);
      }
      else if (!is_border(h, pmesh) && !get(fcm,fh) && !is_border(oh, pmesh) && get(fcm,ofh)){
        set_face(h, boost::graph_traits<PolygonMesh>::null_face(), pmesh);
      }
      else {
        // no face kept
        CGAL_assertion((is_border(h, pmesh) || !get(fcm,fh)) && (is_border(oh, pmesh) || !get(fcm,ofh)));
        // vertices pointing to e must change their halfedge
        if (halfedge(v, pmesh) == oh){
          set_halfedge(v, prev(h, pmesh), pmesh);
        }
        if (halfedge(w, pmesh) == h){
          set_halfedge(w, prev(oh, pmesh), pmesh);
        }
        // shortcut the next pointers as e will be removed
        set_next(prev(h, pmesh), next(oh, pmesh), pmesh);
        set_next(prev(oh, pmesh), next(h, pmesh), pmesh);
        remove_edge(e, pmesh);
      }
    }
    else if (keep_vertex[v]){
      if (halfedge(v, pmesh) == oh){
        set_halfedge(v, prev(h, pmesh), pmesh);
      }
      set_next(prev(h, pmesh), next(oh, pmesh), pmesh);
      remove_edge(e, pmesh);
    }
    else {
      CGAL_assertion(keep_vertex[w]);
      if (halfedge(w, pmesh) == h){
        set_halfedge(w, prev(oh, pmesh), pmesh);
      }
      set_next(prev(oh, pmesh), next(h, pmesh), pmesh);
      remove_edge(e, pmesh);
    }
  }

  face_iterator fb, fe;
  // We now can remove all vertices and faces not marked as kept
  for (boost::tie(fb, fe) = faces(pmesh); fb != fe;){
    face_descriptor f = *fb;
    ++fb;
    if (get(fcm,f) != 1){
      remove_face(f, pmesh);
    }
  }
  vertex_iterator b, e;
  for (boost::tie(b, e) = vertices(pmesh); b != e;){
    vertex_descriptor v = *b;
    ++b;
    if (!keep_vertex[v]){
      remove_vertex(v, pmesh);
    }
  }
}

/*!
* \ingroup keep_connected_components_grp
* keeps the connected components designated by theirs ids in `components_to_keep`,
* and removes the other connected components as well as all isolated vertices.
* The connected component id of a face is given by `fcm`.
*
* \note If the removal of the connected components makes `pmesh` a non-manifold surface,
* then the behavior of this function is undefined.
*
* If `PolygonMesh`
*  has a non writable internal property map
*  for `CGAL::face_index_t` or `CGAL::vertex_index_t` and no `face_index_map` (respectively `vertex_index_map`) is given
*  as a named parameter, then the internal one(s) must be initialized. Otherwise, it will be.
*
* \tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
* \tparam ComponentRange a range of ids convertible to `std::size`
* \tparam FaceComponentMap a model of `ReadWritePropertyMap` with
*         `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
*         `boost::graph_traits<PolygonMesh>::%faces_size_type` as value type.
*
* \param components_to_keep the range of ids of connected components to keep
* \param pmesh the polygon mesh
* \param fcm the property map with indices of components associated to faces in `pmesh`.
*        After calling this function, the values of `fcm` are undefined.
* \param np optional \ref pmp_namedparameters "Named Parameters" described below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap
        , typename NamedParameters>
void keep_connected_components(PolygonMesh& pmesh
                              , const ComponentRange& components_to_keep
                              , const FaceComponentMap& fcm
                              , const NamedParameters& np)
{
  keep_or_remove_connected_components(pmesh, components_to_keep, fcm, true, np);
}

/*!
* \ingroup keep_connected_components_grp
* Removes in `pmesh` the connected components designated by theirs ids
* in `components_to_remove` as well as all isolated vertices.
* The connected component id of a face is given by `fcm`.
*
* \note If the removal of the connected components makes `pmesh` a non-manifold surface,
* then the behavior of this function is undefined.
*
* If `PolygonMesh`
*  has a non writable internal property map
*  for `CGAL::face_index_t` or `CGAL::vertex_index_t` and no `face_index_map` (respectively `vertex_index_map`) is given
*  as a named parameter, then the internal one(s) must be initialized. Otherwise, it will be.
*
*
* \tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
* \tparam ComponentRange a range of ids convertible to `std::size`
* \tparam FaceComponentMap a model of `ReadWritePropertyMap` with
*         `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and
*         `boost::graph_traits<PolygonMesh>::%faces_size_type` as value type.
*
* \param components_to_remove the range of ids of connected components to remove
* \param pmesh the polygon mesh
* \param fcm the property map with indices of components associated to faces in `pmesh`.
*        After calling this function, the values of `fcm` are undefined.
* \param np optional \ref pmp_namedparameters "Named Parameters" described below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap
        , typename NamedParameters>
void remove_connected_components(PolygonMesh& pmesh
                                , const ComponentRange& components_to_remove
                                , const FaceComponentMap& fcm
                                , const NamedParameters& np)
{
  if (components_to_remove.empty()) return;
  keep_or_remove_connected_components(pmesh, components_to_remove, fcm, false, np);
}

/*!
* \ingroup keep_connected_components_grp
*  keeps the connected components not designated by the faces in `components_to_remove`,
*  and removes the other connected components and all isolated vertices.
*
* If `PolygonMesh`
*  has a non writable internal property map
*  for `CGAL::face_index_t` or `CGAL::vertex_index_t` and no `face_index_map` (respectively `vertex_index_map`) is given
*  as a named parameter, then the internal one(s) must be initialized. Otherwise, it will be.
*
* \note If the removal of the connected components makes `pmesh` a non-manifold surface,
* then the behavior of this function is undefined.
*
* \tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
* \tparam FaceRange a range of `boost::graph_traits<PolygonMesh>::%face_descriptor`
*         indicating the connected components to be removed.
*
* \param components_to_remove a face range, including one face or more on each component to be removed
* \param pmesh the polygon mesh
* \param np optional \ref pmp_namedparameters "Named Parameters", amongst those described below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{edge_is_constrained_map} a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
*    \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
*    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
        , typename FaceRange
        , typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void remove_connected_components(PolygonMesh& pmesh
                                , const FaceRange& components_to_remove
                                , const CGAL_PMP_NP_CLASS& np)
{
  if (components_to_remove.empty()) return;
  typedef PolygonMesh PM;
  typedef typename boost::graph_traits<PM>::face_descriptor face_descriptor;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  //FaceIndexMap
  typedef typename Default_face_index_map<CGAL_PMP_NP_CLASS, PolygonMesh>::type FaceIndexMap;
  FaceIndexMap fim =
      CGAL::get_initialized_face_index_map(pmesh, np);

  //vector_property_map
  boost::vector_property_map<std::size_t, FaceIndexMap> face_cc(fim);

  connected_components(pmesh, face_cc, np);

  std::vector<std::size_t> cc_to_remove;
  for(face_descriptor f : components_to_remove)
    cc_to_remove.push_back( face_cc[f] );

  remove_connected_components(pmesh, cc_to_remove, face_cc, np);
}

/*!
* \ingroup keep_connected_components_grp
*  keeps the connected components designated by the faces in `components_to_keep`,
*  and removes the other connected components and all isolated vertices.
*
* If `PolygonMesh`
*  has a non writable internal property map
*  for `CGAL::face_index_t` or `CGAL::vertex_index_t` and no `face_index_map` (respectively `vertex_index_map`) is given
*  as a named parameter, then the internal one(s) must be initialized. Otherwise, it will be.
*
* \note If the removal of the connected components makes `pmesh` a non-manifold surface,
* then the behavior of this function is undefined.
*
* \tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
* \tparam FaceRange a range of `boost::graph_traits<PolygonMesh>::%face_descriptor`
*         indicating the connected components to be kept.
*
* \param pmesh the polygon mesh
* \param components_to_keep a face range, including one face or more on each component to be kept
* \param np optional \ref pmp_namedparameters "Named Parameters", amongst those described below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{edge_is_constrained_map} a property map containing the constrained-or-not status of each edge of `pmesh` \cgalParamEnd
*    \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh` \cgalParamEnd
*    \cgalParamBegin{vertex_index_map} a property map containing the index of each vertex of `pmesh` \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template <typename PolygonMesh
        , typename FaceRange
        , typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void keep_connected_components(PolygonMesh& pmesh
                             , const FaceRange& components_to_keep
                             , const CGAL_PMP_NP_CLASS& np)
{
  typedef PolygonMesh PM;
  typedef typename boost::graph_traits<PM>::face_descriptor face_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  //FaceIndexMap
  typedef typename Default_face_index_map<CGAL_PMP_NP_CLASS, PolygonMesh>::type FaceIndexMap;
  FaceIndexMap fim =
      CGAL::get_initialized_face_index_map(pmesh, np);

  //vector_property_map
  boost::vector_property_map<std::size_t, FaceIndexMap> face_cc(fim);

  connected_components(pmesh, face_cc, np);

  std::vector<std::size_t> cc_to_keep;
  for(face_descriptor f : components_to_keep)
    cc_to_keep.push_back( face_cc[f] );

  keep_connected_components(pmesh, cc_to_keep, face_cc, np);
}

// non-documented overloads so that named parameters can be omitted

template <typename PolygonMesh, typename FaceRange>
void remove_connected_components(PolygonMesh& pmesh
                                , const FaceRange& components_to_remove)
{
  remove_connected_components(pmesh, components_to_remove,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap>
void keep_connected_components(PolygonMesh& pmesh
                              , const ComponentRange& components_to_keep
                              , const FaceComponentMap& fcm)
{
  keep_connected_components(pmesh, components_to_keep, fcm,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

template <typename PolygonMesh
        , typename ComponentRange
        , typename FaceComponentMap>
void remove_connected_components(PolygonMesh& pmesh
                                , const ComponentRange& components_to_remove
                                , const FaceComponentMap& fcm )
{
    remove_connected_components(pmesh, components_to_remove, fcm,
      CGAL::Polygon_mesh_processing::parameters::all_default());
}

template <typename PolygonMesh, typename FaceRange>
void keep_connected_components(PolygonMesh& pmesh
                             , const FaceRange& components_to_keep)
{
  keep_connected_components(pmesh, components_to_keep,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

} // namespace Polygon_mesh_processing

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_POLYGON_MESH_PROCESSING_CONNECTED_COMPONENTS_H
