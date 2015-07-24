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

#define CGAL_DUMP_REMESHING_STEPS

#include <CGAL/Polygon_mesh_processing/internal/remesh_impl.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL {

namespace Polygon_mesh_processing {

/*!
* \ingroup PkgPolygonMeshProcessing
* @brief remeshes a triangulated region of a polygon mesh.
* implements section 6.5.3 "Incremental remeshing" from the PMP book
*
* @tparam PolygonMesh model of `MutableFaceGraph` that 
*         has an internal property map for `CGAL::vertex_point_t`
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `InputIterator`.
* @tparam NamedParameters a sequence of \ref namedparameters
*
* @param pmesh polygon mesh with patches to be remeshed
* @param faces the range of faces defining one patch to remesh
* @param target_edge_length the edge length that is targetted in the remeshed patch
* @param np optional sequence of \ref namedparameters among the ones listed below

* \cgalNamedParamsBegin
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} TODO...
*  \cgalParamEnd
*  \cgalParamBegin{geom_traits} a geometric traits class instance
*  \cgalParamEnd
*  \cgalParamBegin{edge_is_constrained_map}
*  \cgalParamEnd
*  \cgalParamBegin{protect_constraints} If `true`, the edges set as constrained
*     in `edge_is_constrained_map` (or by default the boundary edges)
*     are not splitted nor collapsed during remeshing.
*     Note that around edges that have their length higher than
*     twice `target_edge_length`, remeshing will fail to provide
*     good quality results.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*
* @sa `split_long_edges`
*
*@todo we suppose `faces` describe only one patch. Handle several patches.
*@todo document `number_of_iterations`
*@todo add possibility to provide a functor that projects to a prescribed surface
*/
template<typename PolygonMesh
       , typename FaceRange
       , typename NamedParameters>
void incremental_triangle_based_remeshing(PolygonMesh& pmesh
                                        , const FaceRange& faces
                                        , const double& target_edge_length
                                        , const NamedParameters& np)
{
  typedef PolygonMesh PM;
  using boost::choose_pmap;
  using boost::get_param;
  using boost::choose_param;

  typedef typename GetGeomTraits<PM, NamedParameters>::type GT;

  typedef typename GetVertexPointMap<PM, NamedParameters>::type VPMap;
  VPMap vpmap = choose_pmap(get_param(np, boost::vertex_point),
                            pmesh,
                            boost::vertex_point);

  typedef typename boost::lookup_named_param_def <
      CGAL::edge_is_constrained_t,
      NamedParameters,
      internal::Border_constraint_pmap<PM, FaceRange>//default
    > ::type ECMap;
  ECMap ecmap
    = choose_param(get_param(np, edge_is_constrained),
                   internal::Border_constraint_pmap<PM, FaceRange>(pmesh, faces));

  bool protect = choose_param(get_param(np, protect_constraints), false);

  typename internal::Incremental_remesher<PM, VPMap, GT>
    remesher(pmesh, vpmap, protect);
  remesher.init_faces_remeshing(faces, ecmap);

  unsigned int nb_iterations = choose_param(get_param(np, number_of_iterations), 1);

  double low = 4. / 5. * target_edge_length;
  double high = 4. / 3. * target_edge_length;

#ifdef CGAL_PMP_REMESHING_VERBOSE
  std::cout << std::endl;
  std::cout << "Remeshing (size = " << target_edge_length;
  std::cout << ", #iter = " << nb_iterations << ")..." << std::endl;
#endif

  for (unsigned int i = 0; i < nb_iterations; ++i)
  {
#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout << " * Iteration " << (i + 1) << " *" << std::endl;
#endif

    remesher.split_long_edges(high);
    remesher.collapse_short_edges(low, high);
    remesher.equalize_valences();
    remesher.tangential_relaxation();
    remesher.project_to_surface();

#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout << std::endl;
#endif
  }

#ifdef CGAL_PMP_REMESHING_VERBOSE
  std::cout << "Remeshing done (size = " << target_edge_length;
  std::cout << ", #iter = " << nb_iterations << ")." << std::endl;
#endif
}

template<typename PolygonMesh
       , typename FaceRange>
void incremental_triangle_based_remeshing(PolygonMesh& pmesh
  , const FaceRange& faces
  , const double& target_edge_length)
{
  incremental_triangle_based_remeshing(pmesh,
    faces,
    target_edge_length,
    parameters::all_default());
}

/*!
* \ingroup PkgPolygonMeshProcessing
* @brief splits the edges listed in `edges`
*
* @tparam PolygonMesh model of `MutableFaceGraph` that
*         has an internal property map for `CGAL::vertex_point_t`
*/
template<typename PolygonMesh
       , typename EdgeRange
       , typename NamedParameters>
void split_long_edges(PolygonMesh& pmesh
                    , EdgeRange& edge_range
                    , const double& max_length
                    , const NamedParameters& np)
{
  typedef PolygonMesh PM;
  using boost::choose_pmap;
  using boost::get_param;

  typedef typename GetGeomTraits<PM, NamedParameters>::type GT;
  typedef typename GetVertexPointMap<PM, NamedParameters>::type VPMap;
  VPMap vpmap = choose_pmap(get_param(np, boost::vertex_point),
                            pmesh,
                            boost::vertex_point);

  typename internal::Incremental_remesher<PM, VPMap, GT>
    remesher(pmesh, vpmap, false/*protect constraints*/);

  remesher.split_long_edges(edge_range, max_length, Emptyset_iterator());
}

template<typename PolygonMesh, typename EdgeRange>
void split_long_edges(PolygonMesh& pmesh
                    , EdgeRange& edges
                    , const double& max_length)
{
  split_long_edges(pmesh,
    edges,
    max_length,
    parameters::all_default());
}

//used in the Polyhedron demo
template<typename PolygonMesh
       , typename EdgeRange
       , typename OutputIterator
       , typename NamedParameters>
void split_long_edges(PolygonMesh& pmesh
        , EdgeRange& edge_range
        , const double& max_length
        , OutputIterator out//edges after splitting, all shorter than target_length
        , const NamedParameters& np)
{
  typedef PolygonMesh PM;
  using boost::choose_pmap;
  using boost::get_param;

  typedef typename GetGeomTraits<PM, NamedParameters>::type GT;
  typedef typename GetVertexPointMap<PM, NamedParameters>::type VPMap;
  VPMap vpmap = choose_pmap(get_param(np, boost::vertex_point),
                            pmesh,
                            boost::vertex_point);

  typename internal::Incremental_remesher<PM, VPMap, GT>
    remesher(pmesh, vpmap, false/*protect constraints*/);

  remesher.split_long_edges(edge_range, max_length, out);
}

} //end namespace Polygon_mesh_processing
} //end namespace CGAL

#endif
