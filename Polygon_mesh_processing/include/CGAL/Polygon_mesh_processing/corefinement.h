// Copyright (c) 2016 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_COREFINEMENT_H
#define CGAL_POLYGON_MESH_PROCESSING_COREFINEMENT_H


namespace CGAL {
namespace Polygon_mesh_processing {

// TODO:
//  * parler open/close et si pas closed une seule CC
//  * orientation differente composnate connexes
//  * intersection edge is always on the border of a patch
//  * orientation des polyhedrons
//  * document all_default

/** \ingroup PMP_corefinement_grp
 *
 * indicates if `tm` bounds a volume.
 * See \ref coref_def_subsec for details.
 *
 * @tparam TriangleMesh a model of `FaceGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tm a triangulated surface mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map}
        a property map with the points associated to the vertices of `tm`. The
        point type must be a point from a \cgal Kernel.
     \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \todo in the implementation degenerated faces should be skipt
 */
template <class TriangleMesh, class NamedParameters>
bool does_bound_a_volume(TriangleMesh& tm, const NamedParameters& np);

/**
  * \ingroup PMP_corefinement_grp
  * puts in `O` a triangulated surface mesh \link coref_def_subsec bounding \endlink the union of the volumes
  * bounded by `P` and `Q`.
  * If `O` is one of the input surface meshes, it will be updated to
  * contain the output (in-place operation), otherwise the result will
  * be inserted into `O` without clearing it first.
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(P)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersec(Q)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_bound_a_volume() `CGAL::Polygon_mesh_processing::does_bound_a_volume(P)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_bound_a_volume() `CGAL::Polygon_mesh_processing::does_bound_a_volume(Q)` \endlink
  *
  * @tparam TriangleMesh a model of `FaceGraph`
  * @tparam NamedParametersP a sequence of \ref namedparameters
  * @tparam NamedParametersQ a sequence of \ref namedparameters
  * @tparam NamedParametersO a sequence of \ref namedparameters
  *
  * @param P first input triangulated surface mesh
  * @param Q second input triangulated surface mesh
  * @param O output surface mesh
  * @param np_for_P optional sequence of \ref namedparameters among the ones listed below
  * @param np_for_Q optional sequence of \ref namedparameters among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamBegin{vertex_point_map} a property map with the points associated to the vertices of `P` (`Q`) \cgalParamEnd
  *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
  *     constrained-or-not status of each edge of `P` (`Q`).
  *   \cgalParamBegin{face_index_map} a property map containing the index of each face of `P` (`Q`) \cgalParamEnd
  * \cgalNamedParamsEnd
  *
  * @param np_for_O optional sequence of \ref namedparameters among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamBegin{vertex_point_map} a property map with the points associated to the vertices of `O` \cgalParamEnd
  *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
  *     constrained-or-not status of each edge of `O`. An edge of `O` is constrained
  *     if it is on the intersection of `P` and `Q`, or if the edge corresponds to a
  *     constrained edge in `P` or `Q`.
  * \cgalNamedParamsEnd
  *
  * @return `true` if the output surface mesh is manifold and was put into `O`.
  *         If `false` is returned and `O` was one of the input surface meshes,
  *         then `O` will bound the same volume as the input was but the surface
  *         mesh will be nonetheless \link coref_def_subsec corefined\endlink.
  */
template <class TriangleMesh,
          class NamedParametersP,
          class NamedParametersQ,
          class NamedParametersO,
          >
bool
join(const TriangleMesh& P,
     const TriangleMesh& Q,
           TriangleMesh& O,
     const NamedParametersP& np_for_P,
     const NamedParametersQ& np_for_Q,
     const NamedParametersO& np_for_O);

/**
  * \ingroup PMP_corefinement_grp
  * puts in `O` a triangulated surface mesh \link coref_def_subsec bounding \endlink
  * the intersection of the volumes bounded by `P` and `Q`. Refer to the documentation
  * of \link CGAL::Polygon_mesh_processing::join() `join()` \endlink for
  * details on the parameter.
  */
template <class TriangleMesh,
          class NamedParametersP,
          class NamedParametersQ,
          class NamedParametersO,
          >
bool
intersection(const TriangleMesh& P,
             const TriangleMesh& Q,
                   TriangleMesh& O,
             const NamedParametersP& np_for_P,
             const NamedParametersQ& np_for_Q,
             const NamedParametersO& np_for_O);
/**
  * \ingroup PMP_corefinement_grp
  * puts in `O` a triangulated surface mesh \link coref_def_subsec bounding \endlink
  * the volume bounded by `P` minus the volume bounded by `Q`.
  * Refer to the documentation
  * of \link CGAL::Polygon_mesh_processing::join() `join()` \endlink for
  * details on the parameter.
  */
template <class TriangleMesh,
          class NamedParametersP,
          class NamedParametersQ,
          class NamedParametersO,
          >
bool
difference(const TriangleMesh& P,
           const TriangleMesh& Q,
                 TriangleMesh& O,
           const NamedParametersP& np_for_P,
           const NamedParametersQ& np_for_Q,
           const NamedParametersO& np_for_O);

/**
 * \ingroup PMP_corefinement_grp
 * \link coref_def_subsec corefines \endlink `P` and `Q`. For each input
 * triangulated surface mesh, if a constrained edge is provided, intersection
 * edges will be marked as constrained. If an edge that was marked as
 * constrained is split, its sub-edges will be marked as constrained as well.
 *
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(P)` \endlink
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(Q)` \endlink
 *
 * @tparam TriangleMesh a model of `FaceGraph`
 * @tparam NamedParametersP a sequence of \ref namedparameters
 * @tparam NamedParametersQ a sequence of \ref namedparameters
 *
 * @param P first input triangulated surface mesh
 * @param Q second input triangulated surface mesh
 * @param np_for_P optional sequence of \ref namedparameters among the ones listed below
 * @param np_for_Q optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map} a property map with the points associated to the vertices of `P` (resp. `Q`) \cgalParamEnd
 *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
 *     constrained-or-not status of each edge of `P` (resp. `Q`)\cgalParamEnd
 *   \cgalParamBegin{face_index_map} a property map containing the index of each face of `P` (resp. `Q`) \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
 void
 corefine(const TriangleMesh& P,
          const TriangleMesh& Q,
          const NamedParametersP& np_for_P,
          const NamedParametersQ& np_for_Q);

} }  // end of namespace CGAL::Polygon_mesh_processing

#endif CGAL_POLYGON_MESH_PROCESSING_COREFINEMENT_H
