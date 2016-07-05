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

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Visitor.h>
#include <CGAL/iterator.h>

namespace CGAL {
namespace Polygon_mesh_processing {

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
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm2)` \endlink
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
  * @return `true` if the output surface mesh is manifold and is put into `tm_out`.
  *         If `false` is returned and if `tm_out` is one of the input surface meshes,
  *         then `tm_out` is nonetheless \link coref_def_subsec corefined\endlink
  *         while its bounded volume stays the same as the input one.  */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2,
          class NamedParametersOut>
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
          class NamedParametersOut>
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
          class NamedParametersOut>
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
 * @param throw_on_self_intersection if `true`, for each input triangle mesh,
 *        the set of triangles closed to the intersection of `tm1` and `tm2` will be
 *        checked for self-intersection and `CGAL::Corefinement::Self_intersection_exception`
 *        will be thrown if at least one is found.
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map} a property map with the points associated to the vertices of `tm1` (resp. `tm2`) \cgalParamEnd
 *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
 *     constrained-or-not status of each edge of `tm1` (resp. `tm2`)\cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
 template <class TriangleMesh,
           class NamedParameters1,
           class NamedParameters2>
 void
 corefine(TriangleMesh& tm1,
          TriangleMesh& tm2,
          const NamedParameters1& np1,
          const NamedParameters2& np2,
          const bool throw_on_self_intersection = false)
{
// Vertex point maps
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters1>::const_type Vpm;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters2>::const_type Vpm2;
  CGAL_assertion_code(static const bool same_vpm = )
    boost::is_same<Vpm,Vpm2>::value;
  CGAL_static_assertion(same_vpm);

  Vpm vpm1 = choose_pmap(get_param(np1, boost::vertex_point),
                         tm1,
                         vertex_point);
  Vpm vpm2 = choose_pmap(get_param(np2, boost::vertex_point),
                         tm2,
                         vertex_point);
// Edge is-constrained maps
  typedef typename boost::lookup_named_param_def <
    CGAL::edge_is_constrained_t,
    NamedParameters1,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm1;

  typedef typename boost::lookup_named_param_def <
    CGAL::edge_is_constrained_t,
    NamedParameters1,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm2;

  Ecm1 ecm1 = choose_param( get_param(np1, edge_is_constrained),
                            Corefinement::No_mark<TriangleMesh>() );
  Ecm2 ecm2 = choose_param( get_param(np2, edge_is_constrained),
                            Corefinement::No_mark<TriangleMesh>() );

  typedef Corefinement::Ecm_bind<TriangleMesh, Ecm1, Ecm2> Ecm;

// surface intersection algorithm call
  typedef Corefinement::Default_node_visitor<TriangleMesh> Dnv;
  typedef Corefinement::Default_face_visitor<TriangleMesh> Dfv;
  typedef Corefinement::No_extra_output_from_corefinement<TriangleMesh> Ob;
  typedef Corefinement::Visitor<TriangleMesh,Vpm,Ob,Ecm> Visitor;
  Dnv dnv;
  Dfv dfv;
  Ob ob;
  Ecm ecm(tm1,tm2,ecm1,ecm2);
  Corefinement::Intersection_of_triangle_meshes<TriangleMesh,Vpm,Visitor >
    functor(tm1, tm2, vpm1, vpm2, Visitor(dnv,dfv,ob,ecm));
  functor(CGAL::Emptyset_iterator(), throw_on_self_intersection, true);
}

 template <class TriangleMesh, class NamedParameters1>
 void
 corefine(TriangleMesh& tm1,
          TriangleMesh& tm2,
          const NamedParameters1& np1,
          const bool throw_on_self_intersection = false)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  corefine(tm1, tm2, np1, all_default(), throw_on_self_intersection);
}

 template <class TriangleMesh>
 void
 corefine(TriangleMesh& tm1,
          TriangleMesh& tm2,
          const bool throw_on_self_intersection = false)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  corefine(tm1, tm2, all_default(), all_default(), throw_on_self_intersection);
}

} }  // end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_COREFINEMENT_H
