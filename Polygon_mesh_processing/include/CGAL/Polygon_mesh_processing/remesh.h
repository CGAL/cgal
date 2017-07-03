// Copyright (c) 2015 GeometryFactory (France).
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
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_REMESH_H
#define CGAL_POLYGON_MESH_PROCESSING_REMESH_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>


#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/remesh_impl.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#ifdef CGAL_PMP_REMESHING_VERBOSE
#include <CGAL/Timer.h>
#endif

namespace CGAL {

namespace Polygon_mesh_processing {

/*!
* \ingroup PMP_meshing_grp
* @brief remeshes a triangulated region of a polygon mesh.
* This operation sequentially performs edge splits, edge collapses,
* edge flips, tangential relaxation and projection to the initial surface
* to generate a smooth mesh with a prescribed edge length.
*
* @tparam PolygonMesh model of `MutableFaceGraph`.
*         The descriptor types `boost::graph_traits<PolygonMesh>::%face_descriptor`
*         and `boost::graph_traits<PolygonMesh>::%halfedge_descriptor` must be
*         models of `Hashable`.
*         If `PolygonMesh` has an internal property map for `CGAL::face_index_t`,
*         and no `face_index_map` is given
*         as a named parameter, then the internal one must be initialized
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref namedparameters
*
* @param pmesh a polygon mesh with triangulated surface patches to be remeshed
* @param faces the range of triangular faces defining one or several surface patches to be remeshed
* @param target_edge_length the edge length that is targeted in the remeshed patch.
*        If `0` is passed then only the edge-flip, tangential relaxation, and projection steps will be done.
* @param np optional sequence of \ref namedparameters among the ones listed below
*
* @pre if constraints protection is activated, the constrained edges should
* not be longer than 4/3*`target_edge_length`
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of atomic operations performed (listed in the above description)
*  \cgalParamEnd
*  \cgalParamBegin{edge_is_constrained_map} a property map containing the
*    constrained-or-not status of each edge of `pmesh`. A constrained edge can be split
*    or collapsed, but not flipped, nor its endpoints moved by smoothing.
*    Note that patch boundary edges (i.e. incident to only one face in the range)
*    are always considered as constrained edges.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during remeshing
*  \cgalParamEnd
*  \cgalParamBegin{protect_constraints} If `true`, the edges set as constrained
*     in `edge_is_constrained_map` (or by default the boundary edges)
*     are not split nor collapsed during remeshing.
*     Note that around constrained edges that have their length higher than
*     twice `target_edge_length`, remeshing will fail to provide
*     good quality results. It can even fail to terminate because of cascading vertex
*     insertions.
*  \cgalParamEnd
*  \cgalParamBegin{face_patch_map} a property map with the patch id's associated to the
     faces of `faces`. Instance of a class model of `ReadWritePropertyMap`. It gets
     updated during the remeshing process while new faces are created.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_relaxation_steps} the number of iterations of tangential
*    relaxation that are performed at each iteration of the remeshing process
*  \cgalParamEnd
*  \cgalParamBegin{relax_constraints} If `true`, the end vertices of the edges set as
*    constrained in `edge_is_constrained_map` and boundary edges move along the
*    constrained polylines they belong to.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*
* @sa `split_long_edges()`
*
*@todo add possibility to provide a functor that projects to a prescribed surface
*@todo Deal with exact constructions Kernel. The only thing that makes sense is to
*      guarantee that the output vertices are exactly on the input surface.
*      To do so, we can do every construction in `double`, and use an exact process for
*      projection. For each vertex, the `AABB_tree` would be used in an inexact manner
*      to find the triangle on which projection has to be done. Then, use
*      `CGAL::intersection(triangle, line)` in the exact constructions kernel to
*      get a point which is exactly on the surface.
*
*/
template<typename PolygonMesh
       , typename FaceRange
       , typename NamedParameters>
void isotropic_remeshing(const FaceRange& faces
                       , const double& target_edge_length
                       , PolygonMesh& pmesh
                       , const NamedParameters& np)
{
  if (boost::begin(faces)==boost::end(faces))
    return;

  typedef PolygonMesh PM;
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  using boost::get_param;
  using boost::choose_param;

#ifdef CGAL_PMP_REMESHING_VERBOSE
  std::cout << std::endl;
  CGAL::Timer t;
  std::cout << "Remeshing parameters...";
  std::cout.flush();
  t.start();
#endif

  typedef typename GetGeomTraits<PM, NamedParameters>::type GT;

  typedef typename GetVertexPointMap<PM, NamedParameters>::type VPMap;
  VPMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                             get_property_map(vertex_point, pmesh));

  typedef typename GetFaceIndexMap<PM, NamedParameters>::type FIMap;
  FIMap fimap = choose_param(get_param(np, internal_np::face_index),
                           get_property_map(face_index, pmesh));

  typedef typename boost::lookup_named_param_def <
      internal_np::edge_is_constrained_t,
      NamedParameters,
      internal::Border_constraint_pmap<PM, FaceRange, FIMap>//default
    > ::type ECMap;
  ECMap ecmap = (boost::is_same<ECMap, internal::Border_constraint_pmap<PM, FaceRange, FIMap> >::value)
     //avoid constructing the Border_constraint_pmap if it's not used
    ? choose_param(get_param(np, internal_np::edge_is_constrained)
                 , internal::Border_constraint_pmap<PM, FaceRange, FIMap>(pmesh, faces, fimap))
    : choose_param(get_param(np, internal_np::edge_is_constrained)
                 , internal::Border_constraint_pmap<PM, FaceRange, FIMap>());

  typedef typename boost::lookup_named_param_def <
      internal_np::vertex_is_constrained_t,
      NamedParameters,
      internal::No_constraint_pmap<vertex_descriptor>//default
    > ::type VCMap;
  VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                             internal::No_constraint_pmap<vertex_descriptor>());

  typedef typename boost::lookup_named_param_def <
      internal_np::face_patch_t,
      NamedParameters,
      internal::Connected_components_pmap<PM, ECMap, FIMap>//default
    > ::type FPMap;
  FPMap fpmap = choose_param(
    get_param(np, internal_np::face_patch),
    internal::Connected_components_pmap<PM, ECMap, FIMap>(pmesh, ecmap, fimap,
        boost::is_default_param(get_param(np, internal_np::face_patch))));

  double low = 4. / 5. * target_edge_length;
  double high = 4. / 3. * target_edge_length;

  bool protect = choose_param(get_param(np, internal_np::protect_constraints), false);
  if(protect)
  {
    std::string msg("Isotropic remeshing : protect_constraints cannot be set to");
    msg.append(" true with constraints larger than 4/3 * target_edge_length.");
    msg.append(" Remeshing aborted.");
    CGAL_precondition_msg(
      internal::constraints_are_short_enough(pmesh, ecmap, vpmap, high),
      msg.c_str());
  }

#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << "\rRemeshing parameters done ("<< t.time() <<" sec)" << std::endl;
  std::cout << "Remesher construction...";
  std::cout.flush();
  t.reset(); t.start();
#endif

  typename internal::Incremental_remesher<PM, VPMap, GT, ECMap, VCMap, FPMap, FIMap>
    remesher(pmesh, vpmap, protect, ecmap, vcmap, fpmap, fimap);
  remesher.init_remeshing(faces);

#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
#endif

  unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);
  bool smoothing_1d = choose_param(get_param(np, internal_np::relax_constraints), false);
  unsigned int nb_laplacian = choose_param(get_param(np, internal_np::number_of_relaxation_steps), 1);

#ifdef CGAL_PMP_REMESHING_VERBOSE
  std::cout << std::endl;
  std::cout << "Remeshing (size = " << target_edge_length;
  std::cout << ", #iter = " << nb_iterations << ")..." << std::endl;
  t.reset(); t.start();
#endif

  for (unsigned int i = 0; i < nb_iterations; ++i)
  {
#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout << " * Iteration " << (i + 1) << " *" << std::endl;
#endif
    if (target_edge_length>0)
    {
      remesher.split_long_edges(high);
      remesher.collapse_short_edges(low, high);
    }
    remesher.equalize_valences();
    remesher.tangential_relaxation(smoothing_1d, nb_laplacian);
    remesher.project_to_surface();

#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout << std::endl;
#endif
  }

  remesher.update_constraints_property_map();

#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << "Remeshing done (size = " << target_edge_length;
  std::cout << ", #iter = " << nb_iterations;
  std::cout << ", " << t.time() << " sec )." << std::endl;
#endif
}

template<typename PolygonMesh
       , typename FaceRange>
void isotropic_remeshing(
    const FaceRange& faces
  , const double& target_edge_length
  , PolygonMesh& pmesh)
{
  isotropic_remeshing(
    faces,
    target_edge_length,
    pmesh,
    parameters::all_default());
}

/*!
* \ingroup PMP_meshing_grp
* @brief splits the edges listed in `edges` into sub-edges
* that are not longer than the given threshold `max_length`.
*
* Note this function is useful to split constrained edges before
* calling `isotropic_remeshing()` with protection of constraints
* activated (to match the constrained edge length required by the
* remeshing algorithm to be guaranteed to terminate)
*
* @tparam PolygonMesh model of `MutableFaceGraph` that
*         has an internal property map for `CGAL::vertex_point_t`.
* @tparam EdgeRange range of `boost::graph_traits<PolygonMesh>::%edge_descriptor`,
*   model of `Range`. Its iterator type is `InputIterator`.
* @tparam NamedParameters a sequence of \ref namedparameters
*
* @param pmesh a polygon mesh
* @param edges the range of edges to be split if they are longer than given threshold
* @param max_length the edge length above which an edge from `edges` is split
*        into to sub-edges
* @param np optional \ref namedparameters described below

* \cgalNamedParamsBegin
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*
* @sa `isotropic_remeshing()`
*
*/
template<typename PolygonMesh
       , typename EdgeRange
       , typename NamedParameters>
void split_long_edges(const EdgeRange& edges
                    , const double& max_length
                    , PolygonMesh& pmesh
                    , const NamedParameters& np)
{
  typedef PolygonMesh PM;
  typedef typename boost::graph_traits<PM>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  using boost::choose_param;
  using boost::get_param;

  typedef typename GetGeomTraits<PM, NamedParameters>::type GT;
  typedef typename GetVertexPointMap<PM, NamedParameters>::type VPMap;
  VPMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                             get_property_map(vertex_point, pmesh));

  typedef typename GetFaceIndexMap<PM, NamedParameters>::type FIMap;
  FIMap fimap = choose_param(get_param(np, internal_np::face_index),
                             get_property_map(face_index, pmesh));

  typedef typename boost::lookup_named_param_def <
        internal_np::edge_is_constrained_t,
        NamedParameters,
        internal::No_constraint_pmap<edge_descriptor>//default
      > ::type ECMap;
  ECMap ecmap = choose_param(get_param(np, internal_np::edge_is_constrained),
                             internal::No_constraint_pmap<edge_descriptor>());
  
  typename internal::Incremental_remesher<PM, VPMap, GT, ECMap,
    internal::No_constraint_pmap<vertex_descriptor>,
    internal::Connected_components_pmap<PM, ECMap, FIMap>,
    FIMap
  >
    remesher(pmesh, vpmap, false/*protect constraints*/
           , ecmap
           , internal::No_constraint_pmap<vertex_descriptor>()
           , internal::Connected_components_pmap<PM, ECMap, FIMap>(pmesh, ecmap, fimap, false)
           , fimap
           , false/*need aabb_tree*/);

  remesher.split_long_edges(edges, max_length);
}

template<typename PolygonMesh, typename EdgeRange>
void split_long_edges(const EdgeRange& edges
                    , const double& max_length
                    , PolygonMesh& pmesh)
{
  split_long_edges(edges,
    max_length,
    pmesh,
    parameters::all_default());
}

} //end namespace Polygon_mesh_processing
} //end namespace CGAL

#endif
