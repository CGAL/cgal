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
 * @tparam TriangleMesh a model of `FaceListGraph` and `HalfedgeListGraph`
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
bool does_bound_a_volume(const TriangleMesh& tm, const NamedParameters& np);

/**
  * \ingroup PMP_corefinement_grp
  * puts in `tm_out` a triangulated surface mesh \link coref_def_subsec bounding \endlink the union of the volumes
  * bounded by `tm1` and `tm2`.
  * If `tm_out` is one of the input surface meshes, it will be updated to
  * contain the output (in-place operation), otherwise the result will
  * be inserted into `tm_out` without clearing it first.
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm1)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersec(tm2)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_bound_a_volume() `CGAL::Polygon_mesh_processing::does_bound_a_volume(tm1)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_bound_a_volume() `CGAL::Polygon_mesh_processing::does_bound_a_volume(tm2)` \endlink
  *
  * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`
  * @tparam NamedParameters1 a sequence of \ref namedparameters
  * @tparam NamedParameters2 a sequence of \ref namedparameters
  * @tparam NamedParametersOut a sequence of \ref namedparameters
  *
  * @param tm1 first input triangulated surface mesh
  * @param tm2 second input triangulated surface mesh
  * @param tm_out output surface mesh
  * @param np1 optional sequence of \ref namedparameters among the ones listed below
  * @param np2 optional sequence of \ref namedparameters among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamBegin{vertex_point_map} a property map with the points associated to the vertices of `tm1` (`tm2`) \cgalParamEnd
  *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
  *     constrained-or-not status of each edge of `tm1` (`tm2`).
  *   \cgalParamBegin{face_index_map} a property map containing the index of each face of `tm1` (`tm2`) \cgalParamEnd
  * \cgalNamedParamsEnd
  *
  * @param np_out optional sequence of \ref namedparameters among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamBegin{vertex_point_map} a property map with the points associated to the vertices of `tm_out` \cgalParamEnd
  *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
  *     constrained-or-not status of each edge of `tm_out`. An edge of `tm_out` is constrained
  *     if it is on the intersection of `tm1` and `tm2`, or if the edge corresponds to a
  *     constrained edge in `tm1` or `tm2`.
  * \cgalNamedParamsEnd
  *
  * @return `true` if the output surface mesh is manifold and was put into `tm_out`.
  *         If `false` is returned and if `tm_out` was one of the input surface meshes,
  *         then `tm_out` will bound the same volume as the input but the surface
  *         mesh will be nonetheless \link coref_def_subsec corefined\endlink.
  */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2,
          class NamedParametersOut,
          >
bool
join(const TriangleMesh& tm1,
     const TriangleMesh& tm2,
           TriangleMesh& tm_out,
     const NamedParameters1& np1,
     const NamedParameters2& np2,
     const NamedParametersOut& np_out);

/**
  * \ingroup PMP_corefinement_grp
  * puts in `tm_out` a triangulated surface mesh \link coref_def_subsec bounding \endlink
  * the intersection of the volumes bounded by `tm1` and `tm2`.
  * \copydetails CGAL::Polygon_mesh_processing::join()
  */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2,
          class NamedParametersOut,
          >
bool
intersection(const TriangleMesh& tm1,
             const TriangleMesh& tm2,
                   TriangleMesh& tm_out,
             const NamedParameters1& np1,
             const NamedParameters2& np2,
             const NamedParametersOut& np_out);
/**
  * \ingroup PMP_corefinement_grp
  * puts in `tm_out` a triangulated surface mesh \link coref_def_subsec bounding \endlink
  * the volume bounded by `tm1` minus the volume bounded by `tm2`.
  * \copydetails CGAL::Polygon_mesh_processing::join()
  */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2,
          class NamedParametersOut,
          >
bool
difference(const TriangleMesh& tm1,
           const TriangleMesh& tm2,
                 TriangleMesh& tm_out,
           const NamedParameters1& np1,
           const NamedParameters2& np2,
           const NamedParametersOut& np_out);

/**
 * \ingroup PMP_corefinement_grp
 * \link coref_def_subsec corefines \endlink `tm1` and `tm2`. For each input
 * triangulated surface mesh, if a constrained edge is provided, intersection
 * edges will be marked as constrained. If an edge that was marked as
 * constrained is split, its sub-edges will be marked as constrained as well.
 *
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm1)` \endlink
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm2)` \endlink
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`
 * @tparam NamedParameters1 a sequence of \ref namedparameters
 * @tparam NamedParameters2 a sequence of \ref namedparameters
 *
 * @param tm1 first input triangulated surface mesh
 * @param tm2 second input triangulated surface mesh
 * @param np1 optional sequence of \ref namedparameters among the ones listed below
 * @param np2 optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map} a property map with the points associated to the vertices of `tm1` (resp. `tm2`) \cgalParamEnd
 *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
 *     constrained-or-not status of each edge of `tm1` (resp. `tm2`)\cgalParamEnd
 *   \cgalParamBegin{face_index_map} a property map containing the index of each face of `tm1` (resp. `tm2`) \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
 template <class TriangleMesh,
           class NamedParameters1,
           class NamedParameters2
           >
 void
 corefine(TriangleMesh& tm1,
          TriangleMesh& tm2,
          const NamedParameters1& np1,
          const NamedParameters2& np2);

} }  // end of namespace CGAL::Polygon_mesh_processing

#endif CGAL_POLYGON_MESH_PROCESSING_COREFINEMENT_H
