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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_COREFINEMENT_H
#define CGAL_POLYGON_MESH_PROCESSING_COREFINEMENT_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Visitor.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Face_graph_output_builder.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Output_builder_for_autorefinement.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/iterator.h>

namespace CGAL {

#if !defined(CGAL_NO_DEPRECATED_CODE) && !defined(DOXYGEN_RUNNING)
namespace Corefinement {
using Polygon_mesh_processing::Corefinement::Self_intersection_exception;
}
#endif

namespace Polygon_mesh_processing {

namespace internal {

template <class Kernel, class TriangleMesh, class VD, class Fid_map, class Vpm>
bool recursive_does_bound_a_volume(const TriangleMesh& tm,
                                         Vpm& vpm,
                                         Fid_map& fid_map,
                                         const std::vector<VD>& xtrm_vertices,
                                         boost::dynamic_bitset<>& cc_handled,
                                         const std::vector<std::size_t>& face_cc,
                                         std::size_t xtrm_cc_id,
                                         bool is_parent_outward_oriented)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef Side_of_triangle_mesh<TriangleMesh, Kernel, Vpm> Side_of_tm;
// first check that the orientation of the current cc is consistant with its
// parent cc containing it
  bool new_is_parent_outward_oriented = internal::is_outward_oriented(
         xtrm_vertices[xtrm_cc_id], tm, parameters::vertex_point_map(vpm));
  if (new_is_parent_outward_oriented==is_parent_outward_oriented)
    return false;
  cc_handled.set(xtrm_cc_id);

  std::size_t nb_cc = cc_handled.size();

// get all cc that are inside xtrm_cc_id
  std::vector<face_descriptor> cc_faces;
  BOOST_FOREACH(face_descriptor fd, faces(tm))
  {
    if(face_cc[get(fid_map, fd)]==xtrm_cc_id)
      cc_faces.push_back(fd);
  }

  typename Side_of_tm::AABB_tree aabb_tree(cc_faces.begin(), cc_faces.end(),
                                           tm, vpm);
  Side_of_tm side_of_cc(aabb_tree);

  std::vector<std::size_t> cc_inside;
  for(std::size_t id=0; id<nb_cc; ++id)
  {
    if (cc_handled.test(id)) continue;
    if (side_of_cc(get(vpm,xtrm_vertices[id]))==ON_BOUNDED_SIDE)
      cc_inside.push_back(id);
  }

// check whether we need another recursion for cc inside xtrm_cc_id
  if (!cc_inside.empty())
  {
    std::size_t new_xtrm_cc_id = cc_inside.front();
    boost::dynamic_bitset<> new_cc_handled(nb_cc,0);
    new_cc_handled.set();
    new_cc_handled.reset(new_xtrm_cc_id);
    cc_handled.set(new_xtrm_cc_id);

    std::size_t nb_candidates = cc_inside.size();
    for (std::size_t i=1;i<nb_candidates;++i)
    {
      std::size_t candidate = cc_inside[i];
      if(get(vpm,xtrm_vertices[candidate]).z() >
         get(vpm,xtrm_vertices[new_xtrm_cc_id]).z()) new_xtrm_cc_id=candidate;
      new_cc_handled.reset(candidate);
      cc_handled.set(candidate);
    }

    if ( !internal::recursive_does_bound_a_volume<Kernel>(
           tm, vpm, fid_map, xtrm_vertices, new_cc_handled, face_cc,
           new_xtrm_cc_id, new_is_parent_outward_oriented) ) return false;
  }

// now explore remaining cc included in the same cc as xtrm_cc_id
  boost::dynamic_bitset<> cc_not_handled = ~cc_handled;
  std::size_t new_xtrm_cc_id = cc_not_handled.find_first();
  if (new_xtrm_cc_id == cc_not_handled.npos) return true;

  for (std::size_t candidate = cc_not_handled.find_next(new_xtrm_cc_id);
                   candidate < cc_not_handled.npos;
                   candidate = cc_not_handled.find_next(candidate))
  {
     if(get(vpm,xtrm_vertices[candidate]).z() > get(vpm,xtrm_vertices[new_xtrm_cc_id]).z())
        new_xtrm_cc_id = candidate;
  }

  return internal::recursive_does_bound_a_volume<Kernel>(
            tm, vpm, fid_map, xtrm_vertices, cc_handled, face_cc,
            new_xtrm_cc_id, is_parent_outward_oriented);
}

} //end of namespace internal

namespace Corefinement
{
/** \ingroup PMP_corefinement_grp
 *  Default new-face visitor model of `PMPCorefinementVisitor`.
 *  All of its functions have an empty body. This class can be used as a
 *  base class if only some of the functions of the concept require to be
 *  overridden.
 */
template <class TriangleMesh>
struct Default_visitor;

#ifdef DOXYGEN_RUNNING
/** \ingroup PMP_corefinement_grp
 *  Integer identifiers to refer to a particular Boolean operation in the function `corefine_and_compute_boolean_operations()`.
 */
enum Boolean_operation_type {UNION = 0, INTERSECTION=1,
                             TM1_MINUS_TM2=2, TM2_MINUS_TM1=3, NONE };
#endif
}

/** \ingroup PMP_corefinement_grp
 *
 * indicates if `tm` bounds a volume.
 * See \ref coref_def_subsec for details.
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * @param tm a closed triangulated surface mesh
 * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * @pre `CGAL::is_closed(tm)`
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map}
 *     the property map with the points associated to the vertices of `tm`.
 *     If this parameter is omitted, an internal property map for
 *     `CGAL::vertex_point_t` must be available in `TriangleMesh`
 *   \cgalParamEnd
 *   \cgalParamBegin{face_index_map}
 *     a property map containing the index of each face of `tm`.
 *   \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * \see `CGAL::Polygon_mesh_processing::orient_to_bound_a_volume()`
 */
template <class TriangleMesh, class NamedParameters>
bool does_bound_a_volume(const TriangleMesh& tm, const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters>::const_type Vpm;
  typedef typename GetFaceIndexMap<TriangleMesh,
                                   NamedParameters>::const_type Fid_map;
  typedef typename Kernel_traits<
    typename boost::property_traits<Vpm>::value_type >::Kernel Kernel;

  if (!is_closed(tm)) return false;
  if (!is_triangle_mesh(tm)) return false;

  Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                get_const_property_map(boost::vertex_point, tm));

  Fid_map fid_map = boost::choose_param(boost::get_param(np, internal_np::face_index),
                                        get_const_property_map(boost::face_index, tm));

  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));

  // set the connected component id of each face
  std::size_t nb_cc = connected_components(tm,
                                bind_property_maps(fid_map,make_property_map(face_cc)),
                                parameters::face_index_map(fid_map));

  if (nb_cc == 1)
    return true;

  boost::dynamic_bitset<> cc_handled(nb_cc, 0);

  // extract a vertex with max z coordinate for each connected component
  std::vector<vertex_descriptor> xtrm_vertices(nb_cc, GT::null_vertex());
  BOOST_FOREACH(vertex_descriptor vd, vertices(tm))
  {
    std::size_t cc_id = face_cc[get(fid_map, face(halfedge(vd, tm), tm))];
    if (xtrm_vertices[cc_id]==GT::null_vertex())
      xtrm_vertices[cc_id]=vd;
    else
      if (get(vpm, vd).z()>get(vpm,xtrm_vertices[cc_id]).z())
        xtrm_vertices[cc_id]=vd;
  }

  //extract a vertex with max z amongst all components
  std::size_t xtrm_cc_id=0;
  for(std::size_t id=1; id<nb_cc; ++id)
    if (get(vpm, xtrm_vertices[id]).z()>get(vpm,xtrm_vertices[xtrm_cc_id]).z())
      xtrm_cc_id=id;

  bool is_parent_outward_oriented =
    !internal::is_outward_oriented(xtrm_vertices[xtrm_cc_id], tm, np);

  return internal::recursive_does_bound_a_volume<Kernel>(tm, vpm, fid_map,
                                                         xtrm_vertices,
                                                         cc_handled,
                                                         face_cc,
                                                         xtrm_cc_id,
                                                         is_parent_outward_oriented);
}

/// \cond SKIP_IN_MANUAL
template <class TriangleMesh>
bool does_bound_a_volume(const TriangleMesh& tm)
{
  return does_bound_a_volume(tm, parameters::all_default());
}
/// \endcond

#define CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP(I) \
  typedef typename boost::lookup_named_param_def < \
    internal_np::edge_is_constrained_t, \
    NamedParametersOut##I, \
    Corefinement::No_mark<TriangleMesh> \
  > ::type Ecm_out_##I; \
    Ecm_out_##I ecm_out_##I = \
      boost::choose_param( boost::get_param(cpp11::get<I>(nps_out), internal_np::edge_is_constrained),  \
                           Corefinement::No_mark<TriangleMesh>() );


/**
  * \ingroup PMP_corefinement_grp
  * \link coref_def_subsec corefines \endlink `tm1` and `tm2` and for each triangle mesh `tm_out` passed
  * as an optional in `output` different from `boost::none`, the triangulated surface mesh
  * \link coref_def_subsec bounding \endlink  the result of a particular Boolean operation
  * between the volumes bounded by `tm1` and `tm2` will be put in the corresponding triangle mesh.
  * The positions of the meshes in the array `output` are specific to the Boolean operation to compute
  * and `Corefinement::Boolean_operation_type` encodes and describes the ordering. Constructing the default array
  * means that no Boolean operation will be done. Overwriting a default value will trigger the corresponding
  * operation. In such a case, the address to a valid surface mesh must be provided.
  * The optional named parameters for all output meshes are provided as a `tuple` and follow the same
  * order as the array `output`. A call to `corefine_and_compute_boolean_operations()` with optional
  * named parameters passed for output meshes should be done using `make_tuple()` as the types of
  * named parameters are unspecified.
  *
  * If `tm1` and/or `tm2` are part of the output surface meshes, they will be updated to
  * contain the output (in-place operation), in any other case, the corresponding result will
  * be inserted into the mesh without clearing it first.
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm1)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm2)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_bound_a_volume() `CGAL::Polygon_mesh_processing::does_bound_a_volume(tm1)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_bound_a_volume() `CGAL::Polygon_mesh_processing::does_bound_a_volume(tm2)` \endlink
  *
  * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`
  * @tparam NamedParameters1 a sequence of \ref pmp_namedparameters "Named Parameters"
  * @tparam NamedParameters2 a sequence of \ref pmp_namedparameters "Named Parameters"
  * @tparam NamedParametersOut0 a sequence of \ref pmp_namedparameters "Named Parameters" for computing the union of the volumes bounded by `tm1` and `tm2`
  * @tparam NamedParametersOut1 a sequence of \ref pmp_namedparameters "Named Parameters" for computing the intersection of the volumes bounded by `tm1` and `tm2`
  * @tparam NamedParametersOut2 a sequence of \ref pmp_namedparameters "Named Parameters" for computing the difference of the volumes bounded by `tm1` and `tm2`
  * @tparam NamedParametersOut3 a sequence of \ref pmp_namedparameters "Named Parameters" for computing the difference of the volumes bounded by `tm2` and `tm1`
  *
  * @param tm1 first input triangulated surface mesh
  * @param tm2 second input triangulated surface mesh
  * @param output an array of output surface meshes
  * @param np1 optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  * @param np2 optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamBegin{vertex_point_map}
  *     the property map with the points associated to the vertices of `tm1` (`tm2`).
  *     If this parameter is omitted, an internal property map for
  *     `CGAL::vertex_point_t` should be available in `TriangleMesh`
  *   \cgalParamEnd
  *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
  *     constrained-or-not status of each edge of `tm1` (`tm2`).
  *   \cgalParamEnd
  *   \cgalParamBegin{face_index_map} a property map containing the index of each face of `tm1` (`tm2`).
  *     Note that if the property map is writable, the indices of the faces
  *     of `tm1` and `tm2` will be set after the corefinement is done.
  *   \cgalParamEnd
  *   \cgalParamBegin{visitor} a class model of `PMPCorefinementVisitor`
  *                            that is used to track the creation of new faces  (`np1` only)
  *   \cgalParamEnd
  *   \cgalParamBegin{throw_on_self_intersection} if `true`, for each input triangle mesh,
  *      the set of triangles close to the intersection of `tm1` and `tm2` will be
  *      checked for self-intersection and `CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception`
  *      will be thrown if at least one is found (`np1` only).
  *   \cgalParamEnd
  * \cgalNamedParamsEnd
  *
  * @param nps_out tuple of optional sequences of \ref pmp_namedparameters "Named Parameters" each among the ones listed below
  *        (`tm_out` being used to refer to the output surface mesh in `output` corresponding to a given named parameter sequence)
  *
  * \cgalNamedParamsBegin
  *   \cgalParamBegin{vertex_point_map}
  *     the property map with the points associated to the vertices of `tm_out`.
  *     If this parameter is omitted, an internal property map for
  *     `CGAL::vertex_point_t` must be available in `TriangleMesh`
  *   \cgalParamEnd
  *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
  *     constrained-or-not status of each edge of `tm_out`. An edge of `tm_out` is constrained
  *     if it is on the intersection of `tm1` and `tm2`, or if the edge corresponds to a
  *     constrained edge in `tm1` or `tm2`.
  *   \cgalParamEnd
  * \cgalNamedParamsEnd
  *
  * @return an array filled as follows: for each operation computed, the position in the array
  *         will contain `true` iff the output surface mesh is manifold, and it is put in the surface mesh
  *         at the same position as in `output`. Note that if an output surface mesh also was
  *         an input mesh but the output operation was generating a non-manifold mesh, the surface mesh
  *         will only be corefined.
  */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2,
          class NamedParametersOut0,
          class NamedParametersOut1,
          class NamedParametersOut2,
          class NamedParametersOut3>
cpp11::array<bool,4>
corefine_and_compute_boolean_operations(
        TriangleMesh& tm1,
        TriangleMesh& tm2,
  const cpp11::array< boost::optional<TriangleMesh*>,4>& output,
  const NamedParameters1& np1,
  const NamedParameters2& np2,
  const cpp11::tuple<NamedParametersOut0,
                     NamedParametersOut1,
                     NamedParametersOut2,
                     NamedParametersOut3>& nps_out)
{
  const bool throw_on_self_intersection =
    boost::choose_param(boost::get_param(np1, internal_np::throw_on_self_intersection), false);

// Vertex point maps
  //for input meshes
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters1>::type Vpm;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters2>::type Vpm2;
  CGAL_USE_TYPE(Vpm2);
  CGAL_assertion_code(
    static const bool same_vpm = (boost::is_same<Vpm,Vpm2>::value); )
  CGAL_static_assertion(same_vpm);

  Vpm vpm1 = boost::choose_param(boost::get_param(np1, internal_np::vertex_point),
                                 get_property_map(boost::vertex_point, tm1));

  Vpm vpm2 = boost::choose_param(boost::get_param(np2, internal_np::vertex_point),
                                 get_property_map(boost::vertex_point, tm2));

  typedef typename boost::property_traits<Vpm>::value_type Point_3;

  // for output meshes: here we have to use a trick so that if for a specific output
  // that is not requested, the default vpm does not have the same value type as the
  // input map, a dummy default vpm is used so that calls to get/put can be compiled
  // (even if not used).
  typedef cpp11::tuple<
    Corefinement::TweakedGetVertexPointMap<Point_3, NamedParametersOut0, TriangleMesh>,
    Corefinement::TweakedGetVertexPointMap<Point_3, NamedParametersOut1, TriangleMesh>,
    Corefinement::TweakedGetVertexPointMap<Point_3, NamedParametersOut2, TriangleMesh>,
    Corefinement::TweakedGetVertexPointMap<Point_3, NamedParametersOut3, TriangleMesh>
  > Vpm_out_tuple_helper;

  typedef cpp11::tuple<
    boost::optional< typename cpp11::tuple_element<0, Vpm_out_tuple_helper>::type::type >,
    boost::optional< typename cpp11::tuple_element<1, Vpm_out_tuple_helper>::type::type >,
    boost::optional< typename cpp11::tuple_element<2, Vpm_out_tuple_helper>::type::type >,
    boost::optional< typename cpp11::tuple_element<3, Vpm_out_tuple_helper>::type::type >
  > Vpm_out_tuple;

  Vpm_out_tuple vpm_out_tuple(
    Corefinement::get_vpm<Point_3>(cpp11::get<0>(nps_out), output[0],
                                   typename cpp11::tuple_element<0, Vpm_out_tuple_helper>::type::Use_default_tag()),
    Corefinement::get_vpm<Point_3>(cpp11::get<1>(nps_out), output[1],
                                   typename cpp11::tuple_element<1, Vpm_out_tuple_helper>::type::Use_default_tag()),
    Corefinement::get_vpm<Point_3>(cpp11::get<2>(nps_out), output[2],
                                   typename cpp11::tuple_element<2, Vpm_out_tuple_helper>::type::Use_default_tag()),
    Corefinement::get_vpm<Point_3>(cpp11::get<3>(nps_out), output[3],
                                   typename cpp11::tuple_element<3, Vpm_out_tuple_helper>::type::Use_default_tag())
  );

  if (&tm1==&tm2)
  {
    // for now edges in a coplanar patch are not constrained so there is nothing to constrained here
    // \todo marked edges from input to output are not ported

    if (output[Corefinement::UNION] != boost::none)
      if (&tm1 != *output[Corefinement::UNION])
        copy_face_graph(tm1,
            *(*output[Corefinement::UNION]),
                        parameters::vertex_point_map(vpm1),
                        parameters::vertex_point_map(*cpp11::get<Corefinement::UNION>(vpm_out_tuple)));
    if (output[Corefinement::INTERSECTION] != boost::none)
      if (&tm1 != *output[Corefinement::INTERSECTION])
        copy_face_graph(tm1,
                        *(*output[Corefinement::INTERSECTION]),
                        parameters::vertex_point_map(vpm1),
                        parameters::vertex_point_map(*cpp11::get<Corefinement::INTERSECTION>(vpm_out_tuple)));
                        

    if (output[Corefinement::TM1_MINUS_TM2] != boost::none)
      if (&tm1 == *output[Corefinement::TM1_MINUS_TM2])
        clear(tm1);

    if (output[Corefinement::TM2_MINUS_TM1] != boost::none)
      if (&tm1 == *output[Corefinement::TM2_MINUS_TM1])
        clear(tm1);

    return CGAL::make_array(true, true, true, true);
  }

  // handle case of empty meshes (isolated vertices are ignored)
  if (faces(tm1).empty())
  {
    if(faces(tm2).empty())
    {
      for (int i=0; i<4; ++i)
        if (output[i] != boost::none)
          clear(*(*output[i]));
      return CGAL::make_array(true, true, true, true);
    }
    // tm2 is not empty
    if (output[Corefinement::UNION] != boost::none)
      if (&tm2 != *output[Corefinement::UNION])
        copy_face_graph(tm2,
                        *(*output[Corefinement::UNION]),
                        parameters::vertex_point_map(vpm2),
                        parameters::vertex_point_map(*cpp11::get<Corefinement::UNION>(vpm_out_tuple)));
    if (output[Corefinement::INTERSECTION] != boost::none)
      clear(*(*output[Corefinement::INTERSECTION]));
    if (output[Corefinement::TM1_MINUS_TM2] != boost::none)
      clear(*(*output[Corefinement::TM1_MINUS_TM2]));
    if (output[Corefinement::TM2_MINUS_TM1] != boost::none)
      if (&tm2 != *output[Corefinement::TM2_MINUS_TM1])
        copy_face_graph(tm2,
                        *(*output[Corefinement::TM2_MINUS_TM1]),
                        parameters::vertex_point_map(vpm2),
                        parameters::vertex_point_map(*cpp11::get<Corefinement::TM2_MINUS_TM1>(vpm_out_tuple)));
    return CGAL::make_array(true, true, true, true);
  }
  else
    if (faces(tm2).empty())
    {
      // tm1 is not empty
      if (output[Corefinement::UNION] != boost::none)
        if (&tm1 != *output[Corefinement::UNION])
          copy_face_graph(tm1,
                          *(*output[Corefinement::UNION]),
                          parameters::vertex_point_map(vpm1),
                          parameters::vertex_point_map(*cpp11::get<Corefinement::UNION>(vpm_out_tuple)));
      if (output[Corefinement::INTERSECTION] != boost::none)
        clear(*(*output[Corefinement::INTERSECTION]));
      if (output[Corefinement::TM2_MINUS_TM1] != boost::none)
        clear(*(*output[Corefinement::TM2_MINUS_TM1]));
      if (output[Corefinement::TM1_MINUS_TM2] != boost::none)
        if (&tm1 != *output[Corefinement::TM1_MINUS_TM2])
          copy_face_graph(tm1,
                          *(*output[Corefinement::TM1_MINUS_TM2]),
                          parameters::vertex_point_map(vpm1),
                          parameters::vertex_point_map(*cpp11::get<Corefinement::TM1_MINUS_TM2>(vpm_out_tuple)));
      return CGAL::make_array(true, true, true, true);
    }

// Edge is-constrained maps
  //for input meshes
  typedef typename boost::lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters1,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm1;

  typedef typename boost::lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters2,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm2;

  Ecm1 ecm1 = boost::choose_param( boost::get_param(np1, internal_np::edge_is_constrained),
                                   Corefinement::No_mark<TriangleMesh>() );
  Ecm2 ecm2 = boost::choose_param( boost::get_param(np2, internal_np::edge_is_constrained),
                                   Corefinement::No_mark<TriangleMesh>() );

  typedef Corefinement::Ecm_bind<TriangleMesh, Ecm1, Ecm2> Ecm_in;

  //for output meshes
  CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP(0)
  CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP(1)
  CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP(2)
  CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP(3)

  // In the current version all types must be the same so an array would be fine
  typedef cpp11::tuple<Ecm_out_0, Ecm_out_1, Ecm_out_2, Ecm_out_3>
                                                            Edge_mark_map_tuple;

// Face index point maps
  typedef typename GetFaceIndexMap<TriangleMesh,
                                   NamedParameters1>::type Fid_map;
  typedef typename GetFaceIndexMap<TriangleMesh,
                                   NamedParameters2>::type Fid_map2;
  CGAL_USE_TYPE(Fid_map2);
  CGAL_assertion_code(
    static const bool same_fidmap = (boost::is_same<Fid_map,Fid_map2>::value);)
  CGAL_static_assertion(same_fidmap);

  Fid_map fid_map1 = boost::choose_param(boost::get_param(np1, internal_np::face_index),
                                        get_property_map(boost::face_index, tm1));
  Fid_map fid_map2 = boost::choose_param(boost::get_param(np2, internal_np::face_index),
                                         get_property_map(boost::face_index, tm2));
// User visitor
  typedef typename boost::lookup_named_param_def <
    internal_np::graph_visitor_t,
    NamedParameters1,
    Corefinement::Default_visitor<TriangleMesh>//default
  > ::type User_visitor;
  User_visitor uv( boost::choose_param( boost::get_param(np1, internal_np::graph_visitor),
                   Corefinement::Default_visitor<TriangleMesh>() ) );

  // surface intersection algorithm call
  typedef Corefinement::Face_graph_output_builder<TriangleMesh,
                                                  Vpm,
                                                  Vpm_out_tuple,
                                                  Fid_map,
                                                  Default,
                                                  Ecm_in,
                                                  Edge_mark_map_tuple,
                                                  User_visitor> Ob;

  typedef Corefinement::Surface_intersection_visitor_for_corefinement<
    TriangleMesh, Vpm, Ob, Ecm_in, User_visitor> Algo_visitor;
  Ecm_in ecm_in(tm1,tm2,ecm1,ecm2);
  Edge_mark_map_tuple ecms_out(ecm_out_0, ecm_out_1, ecm_out_2, ecm_out_3);
  Ob ob(tm1, tm2, vpm1, vpm2, fid_map1, fid_map2, ecm_in,
        vpm_out_tuple, ecms_out, uv, output);

  Corefinement::Intersection_of_triangle_meshes<TriangleMesh, Vpm, Algo_visitor >
    functor(tm1, tm2, vpm1, vpm2, Algo_visitor(uv,ob,ecm_in));
  functor(CGAL::Emptyset_iterator(), throw_on_self_intersection, true);


  return CGAL::make_array(ob.union_is_valid(),
                          ob.intersection_is_valid(),
                          ob.tm1_minus_tm2_is_valid(),
                          ob.tm2_minus_tm1_is_valid());
}

template <class TriangleMesh>
cpp11::array<bool,4>
corefine_and_compute_boolean_operations(
        TriangleMesh& tm1,
        TriangleMesh& tm2,
  const cpp11::array< boost::optional<TriangleMesh*>,4>& output)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_boolean_operations(tm1, tm2, output,
                                                 all_default(), all_default(),
                                                 cpp11::make_tuple(all_default(), all_default(),
                                                                   all_default(), all_default()));
}

template <class TriangleMesh, class NamedParameters1>
cpp11::array<bool,4>
corefine_and_compute_boolean_operations(
        TriangleMesh& tm1,
        TriangleMesh& tm2,
  const cpp11::array< boost::optional<TriangleMesh*>,4>& output,
  const NamedParameters1& np1)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_boolean_operations(tm1, tm2, output,
                                                 np1, all_default(),
                                                 cpp11::make_tuple(all_default(), all_default(),
                                                                   all_default(), all_default()));
}

template <class TriangleMesh, class NamedParameters1, class NamedParameters2>
cpp11::array<bool,4>
corefine_and_compute_boolean_operations(
        TriangleMesh& tm1,
        TriangleMesh& tm2,
  const cpp11::array< boost::optional<TriangleMesh*>,4>& output,
  const NamedParameters1& np1,
  const NamedParameters2& np2)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_boolean_operations(tm1, tm2, output,
                                                 np1, np2,
                                                 cpp11::make_tuple(all_default(), all_default(),
                                                                   all_default(), all_default()));
}

#undef CGAL_COREF_SET_OUTPUT_VERTEX_POINT_MAP
#undef CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP


/**
  * \ingroup PMP_corefinement_grp
  * \link coref_def_subsec corefines \endlink `tm1` and `tm2` and
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
  * @tparam NamedParameters1 a sequence of \ref pmp_namedparameters "Named Parameters"
  * @tparam NamedParameters2 a sequence of \ref pmp_namedparameters "Named Parameters"
  * @tparam NamedParametersOut a sequence of \ref pmp_namedparameters "Named Parameters"
  *
  * @param tm1 first input triangulated surface mesh
  * @param tm2 second input triangulated surface mesh
  * @param tm_out output surface mesh
  * @param np1 optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  * @param np2 optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamBegin{vertex_point_map}
  *     the property map with the points associated to the vertices of `tm1` (`tm2`).
  *     If this parameter is omitted, an internal property map for
  *     `CGAL::vertex_point_t` must be available in `TriangleMesh`
  *   \cgalParamEnd
  *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
  *     constrained-or-not status of each edge of `tm1` (`tm2`).
  *   \cgalParamEnd
  *   \cgalParamBegin{face_index_map} a property map containing the index of each face of `tm1` (`tm2`).
  *     Note that if the property map is writable, the indices of the faces
  *     of `tm1` and `tm2` will be set after the corefinement is done.
  *   \cgalParamEnd
  *   \cgalParamBegin{visitor} a class model of `PMPCorefinementVisitor`
  *                            that is used to track the creation of new faces  (`np1` only)
  *   \cgalParamEnd
  *   \cgalParamBegin{throw_on_self_intersection} if `true`, for each input triangle mesh,
  *      the set of triangles close to the intersection of `tm1` and `tm2` will be
  *      checked for self-intersection and `CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception`
  *      will be thrown if at least one is found (`np1` only).
  *   \cgalParamEnd
  * \cgalNamedParamsEnd
  *
  * @param np_out optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamBegin{vertex_point_map}
  *     the property map with the points associated to the vertices of `tm_out`.
  *     If this parameter is omitted, an internal property map for
  *     `CGAL::vertex_point_t` must be available in `TriangleMesh`
  *   \cgalParamEnd
  *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
  *     constrained-or-not status of each edge of `tm_out`. An edge of `tm_out` is constrained
  *     if it is on the intersection of `tm1` and `tm2`, or if the edge corresponds to a
  *     constrained edge in `tm1` or `tm2`.
  *   \cgalParamEnd
  * \cgalNamedParamsEnd
  *
  * @return `true` if the output surface mesh is manifold and is put into `tm_out`.
  *         If `false` is returned and if `tm_out` is one of the input surface meshes,
  *         then `tm_out` is only corefined.  */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2,
          class NamedParametersOut>
bool
corefine_and_compute_union(      TriangleMesh& tm1,
                                 TriangleMesh& tm2,
                                 TriangleMesh& tm_out,
                           const NamedParameters1& np1,
                           const NamedParameters2& np2,
                           const NamedParametersOut& np_out)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  cpp11::array< boost::optional<TriangleMesh*>,4> output;
  output[Corefinement::UNION]=&tm_out;

  return
   corefine_and_compute_boolean_operations(tm1, tm2, output, np1, np2,
                                           cpp11::make_tuple(np_out,
                                                             all_default(),
                                                             all_default(),
                                                             all_default()))
                                                                [Corefinement::UNION];
}

/**
  * \ingroup PMP_corefinement_grp
  * \link coref_def_subsec corefines \endlink `tm1` and `tm2` and
  * puts in `tm_out` a triangulated surface mesh \link coref_def_subsec bounding \endlink
  * the intersection of the volumes bounded by `tm1` and `tm2`.
  * \copydetails CGAL::Polygon_mesh_processing::corefine_and_compute_union()
  */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2,
          class NamedParametersOut>
bool
corefine_and_compute_intersection(      TriangleMesh& tm1,
                                        TriangleMesh& tm2,
                                        TriangleMesh& tm_out,
                                  const NamedParameters1& np1,
                                  const NamedParameters2& np2,
                                  const NamedParametersOut& np_out)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  cpp11::array< boost::optional<TriangleMesh*>,4> output;
  output[Corefinement::INTERSECTION]=&tm_out;

  return
    corefine_and_compute_boolean_operations(tm1, tm2, output, np1, np2,
                                            cpp11::make_tuple(all_default(),
                                                              np_out,
                                                              all_default(),
                                                              all_default()))
                                                                [Corefinement::INTERSECTION];
}

/**
  * \ingroup PMP_corefinement_grp
  * \link coref_def_subsec corefines \endlink `tm1` and `tm2` and
  * puts in `tm_out` a triangulated surface mesh \link coref_def_subsec bounding \endlink
  * the volume bounded by `tm1` minus the volume bounded by `tm2`.
  * \copydetails CGAL::Polygon_mesh_processing::corefine_and_compute_union()
  */
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2,
          class NamedParametersOut>
bool
corefine_and_compute_difference(      TriangleMesh& tm1,
                                      TriangleMesh& tm2,
                                      TriangleMesh& tm_out,
                                const NamedParameters1& np1,
                                const NamedParameters2& np2,
                                const NamedParametersOut& np_out)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  using namespace CGAL::Polygon_mesh_processing::Corefinement;
  cpp11::array< boost::optional<TriangleMesh*>,4> output;
  output[TM1_MINUS_TM2]=&tm_out;

  return
    corefine_and_compute_boolean_operations(tm1, tm2, output, np1, np2,
                                            cpp11::make_tuple(all_default(),
                                                              all_default(),
                                                              np_out,
                                                              all_default()))
                                                                [TM1_MINUS_TM2];
}

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
 * @tparam NamedParameters1 a sequence of \ref pmp_namedparameters "Named Parameters"
 * @tparam NamedParameters2 a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * @param tm1 first input triangulated surface mesh
 * @param tm2 second input triangulated surface mesh
 * @param np1 optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 * @param np2 optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map}
 *     the property map with the points associated to the vertices of `tm1` (`tm2`).
 *     If this parameter is omitted, an internal property map for
 *     `CGAL::vertex_point_t` must be available in `TriangleMesh`
 *   \cgalParamEnd
 *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
 *     constrained-or-not status of each edge of `tm1` (`tm2`)
 *   \cgalParamEnd
 *   \cgalParamBegin{visitor} a class model of `PMPCorefinementVisitor`
 *                            that is used to track the creation of new faces (`np1` only)
 *   \cgalParamEnd
 *   \cgalParamBegin{throw_on_self_intersection} if `true`, for each input triangle mesh,
 *      the set of triangles close to the intersection of `tm1` and `tm2` will be
 *      checked for self-intersection and `CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception`
 *      will be thrown if at least one is found (`np1` only).
 *   \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
 template <class TriangleMesh,
           class NamedParameters1,
           class NamedParameters2>
 void
 corefine(      TriangleMesh& tm1,
                TriangleMesh& tm2,
          const NamedParameters1& np1,
          const NamedParameters2& np2)
{
  const bool throw_on_self_intersection =
    boost::choose_param(boost::get_param(np1, internal_np::throw_on_self_intersection), false);

// Vertex point maps
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters1>::type Vpm;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters2>::type Vpm2;
  CGAL_USE_TYPE(Vpm2);
  CGAL_assertion_code(
    static const bool same_vpm = (boost::is_same<Vpm,Vpm2>::value);)
  CGAL_static_assertion(same_vpm);

  Vpm vpm1 = boost::choose_param(boost::get_param(np1, internal_np::vertex_point),
                                 get_property_map(boost::vertex_point, tm1));

  Vpm vpm2 = boost::choose_param(boost::get_param(np2, internal_np::vertex_point),
                                 get_property_map(boost::vertex_point, tm2));

// Edge is-constrained maps
  typedef typename boost::lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters1,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm1;

  typedef typename boost::lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters2,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm2;

  Ecm1 ecm1 = boost::choose_param( boost::get_param(np1, internal_np::edge_is_constrained),
                                   Corefinement::No_mark<TriangleMesh>() );
  Ecm2 ecm2 = boost::choose_param( boost::get_param(np2, internal_np::edge_is_constrained),
                                   Corefinement::No_mark<TriangleMesh>() );

  typedef Corefinement::Ecm_bind<TriangleMesh, Ecm1, Ecm2> Ecm;

  if (&tm1==&tm2)
  {
    Corefinement::mark_all_edges(tm1, ecm1);
    Corefinement::mark_all_edges(tm2, ecm2);
    return;
  }

  // User visitor
  typedef typename boost::lookup_named_param_def <
    internal_np::graph_visitor_t,
    NamedParameters1,
    Corefinement::Default_visitor<TriangleMesh>//default
  > ::type User_visitor;
  User_visitor uv( boost::choose_param( boost::get_param(np1, internal_np::graph_visitor),
                   Corefinement::Default_visitor<TriangleMesh>() ) );

// surface intersection algorithm call
  typedef Corefinement::No_extra_output_from_corefinement<TriangleMesh> Ob;
  typedef Corefinement::Surface_intersection_visitor_for_corefinement<
    TriangleMesh, Vpm, Ob, Ecm, User_visitor> Algo_visitor;
  Ob ob;
  Ecm ecm(tm1,tm2,ecm1,ecm2);
  Corefinement::Intersection_of_triangle_meshes<TriangleMesh, Vpm, Algo_visitor>
    functor(tm1, tm2, vpm1, vpm2, Algo_visitor(uv,ob,ecm));
  functor(CGAL::Emptyset_iterator(), throw_on_self_intersection, true);
}

namespace experimental {
/**
 * \ingroup PMP_corefinement_grp
 * \link coref_def_subsec autorefines \endlink `tm`. Refines a triangle mesh
 * so that no triangles intersects in their interior.
 * Self-intersection edges will be marked as constrained. If an edge that was marked as
 * constrained is split, its sub-edges will be marked as constrained as well.
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tm input triangulated surface mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map}
 *     the property map with the points associated to the vertices of `tm`.
 *     If this parameter is omitted, an internal property map for
 *     `CGAL::vertex_point_t` must be available in `TriangleMesh`
 *   \cgalParamEnd
 *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
 *     constrained-or-not status of each edge of `tm`
 *   \cgalParamEnd
 *   \cgalParamBegin{visitor} a class model of `PMPCorefinementVisitor`
 *                            that is used to track the creation of new faces
 *   \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
 template <class TriangleMesh,
           class NamedParameters>
 void
 autorefine(      TriangleMesh& tm,
            const NamedParameters& np)
{
// Vertex point maps
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters>::type Vpm;

  Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                get_property_map(boost::vertex_point, tm));

// Edge is-constrained maps
  typedef typename boost::lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm;


  Ecm ecm = boost::choose_param( boost::get_param(np, internal_np::edge_is_constrained),
                                 Corefinement::No_mark<TriangleMesh>() );

// User visitor
  typedef typename boost::lookup_named_param_def <
    internal_np::graph_visitor_t,
    NamedParameters,
    Corefinement::Default_visitor<TriangleMesh>//default
  > ::type User_visitor;
  User_visitor uv( boost::choose_param( boost::get_param(np, internal_np::graph_visitor),
                   Corefinement::Default_visitor<TriangleMesh>() ) );


// surface intersection algorithm call
  typedef Corefinement::No_extra_output_from_corefinement<TriangleMesh> Ob;
  typedef Corefinement::Surface_intersection_visitor_for_corefinement<
    TriangleMesh, Vpm, Ob, Ecm, User_visitor,true> Algo_visitor;
  Ob ob;

  Corefinement::Intersection_of_triangle_meshes<TriangleMesh, Vpm, Algo_visitor>
    functor(tm, vpm, Algo_visitor(uv,ob,ecm) );

  functor(CGAL::Emptyset_iterator(), true);
}

/**
 * \ingroup PMP_corefinement_grp
 * Removes self-intersections in `tm` by \link coref_def_subsec autorefining \endlink `tm`,
 * removing extra patches, and stitching self-intersection edges.
 * Self-intersection edges will be marked as constrained. If an edge that was marked as
 * constrained is split, its sub-edges will be marked as constrained as well.
 * \return `true` if all self-intersections were fixed and `false` otherwise.
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tm input triangulated surface mesh
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map}
 *     the property map with the points associated to the vertices of `tm`.
 *     If this parameter is omitted, an internal property map for
 *     `CGAL::vertex_point_t` must be available in `TriangleMesh`
 *   \cgalParamEnd
 *   \cgalParamBegin{edge_is_constrained_map} a property map containing the
 *     constrained-or-not status of each edge of `tm`
 *   \cgalParamEnd
 *   \cgalParamBegin{face_index_map} a property map containing the index of each face of `tm` \cgalParamEnd
 *   \cgalParamBegin{visitor} a class model of `PMPCorefinementVisitor`
 *                            that is used to track the creation of new faces
 *   \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
 template <class TriangleMesh,
           class NamedParameters>
 bool
 autorefine_and_remove_self_intersections(      TriangleMesh& tm,
                                          const NamedParameters& np)
{
// Vertex point maps
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters>::type Vpm;
  Vpm vpm = boost::choose_param(boost::get_param(np, internal_np::vertex_point),
                                get_property_map(boost::vertex_point, tm));
// Face index map
  typedef typename GetFaceIndexMap<TriangleMesh,
                                   NamedParameters>::type Fid_map;
  Fid_map fid_map = boost::choose_param(boost::get_param(np, internal_np::face_index),
                                        get_property_map(boost::face_index, tm));
// Edge is-constrained maps
  typedef typename boost::lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm;
  Ecm ecm = boost::choose_param( boost::get_param(np, internal_np::edge_is_constrained),
                                 Corefinement::No_mark<TriangleMesh>() );
// User visitor
  typedef typename boost::lookup_named_param_def <
    internal_np::graph_visitor_t,
    NamedParameters,
    Corefinement::Default_visitor<TriangleMesh>//default
  > ::type User_visitor;
  User_visitor uv( boost::choose_param( boost::get_param(np, internal_np::graph_visitor),
                   Corefinement::Default_visitor<TriangleMesh>() ) );

// surface intersection algorithm call
  typedef Corefinement::Output_builder_for_autorefinement<TriangleMesh,
                                                          Vpm,
                                                          Fid_map,
                                                          Ecm,
                                                          Default > Ob;

  typedef Corefinement::Surface_intersection_visitor_for_corefinement<
    TriangleMesh, Vpm, Ob, Ecm, User_visitor,true> Algo_visitor;
  Ob ob(tm, vpm, fid_map, ecm);

  Corefinement::Intersection_of_triangle_meshes<TriangleMesh, Vpm, Algo_visitor>
    functor(tm, vpm, Algo_visitor(uv,ob,ecm) );

  functor(CGAL::Emptyset_iterator(), true);

  return ob.all_self_intersection_fixed();
}

}// end of namespace experimental

// overload with default named parameters
///// corefine_and_compute_union /////
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool
corefine_and_compute_union(      TriangleMesh& tm1,
                                 TriangleMesh& tm2,
                                 TriangleMesh& tm_out,
                           const NamedParameters1& np1,
                           const NamedParameters2& np2)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_union(tm1, tm2, tm_out,
                                    np1, np2, all_default());
}

template <class TriangleMesh,
          class NamedParameters1>
bool
corefine_and_compute_union(      TriangleMesh& tm1,
                                 TriangleMesh& tm2,
                                 TriangleMesh& tm_out,
                           const NamedParameters1& np1)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_union(tm1, tm2, tm_out,
                                     np1, all_default(), all_default());
}

template <class TriangleMesh>
bool
corefine_and_compute_union(TriangleMesh& tm1,
                           TriangleMesh& tm2,
                           TriangleMesh& tm_out)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_union(tm1, tm2, tm_out,
                                    all_default(), all_default(), all_default());
}

///// corefine_and_compute_intersection /////
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool
corefine_and_compute_intersection(       TriangleMesh& tm1,
                                         TriangleMesh& tm2,
                                         TriangleMesh& tm_out,
                                  const  NamedParameters1& np1,
                                  const  NamedParameters2& np2)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_intersection(tm1, tm2, tm_out,
                                           np1, np2, all_default());
}

template <class TriangleMesh,
          class NamedParameters1>
bool
corefine_and_compute_intersection(      TriangleMesh& tm1,
                                        TriangleMesh& tm2,
                                        TriangleMesh& tm_out,
                                  const NamedParameters1& np1)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_intersection(tm1, tm2, tm_out,
                                           np1, all_default(), all_default());
}

template <class TriangleMesh>
bool
corefine_and_compute_intersection(TriangleMesh& tm1,
                                  TriangleMesh& tm2,
                                  TriangleMesh& tm_out)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_intersection(tm1, tm2, tm_out,
                                           all_default(), all_default(), all_default());
}

///// difference /////
template <class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2>
bool
corefine_and_compute_difference(      TriangleMesh& tm1,
                                      TriangleMesh& tm2,
                                      TriangleMesh& tm_out,
                                const NamedParameters1& np1,
                                const NamedParameters2& np2)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_difference(tm1, tm2, tm_out,
                                         np1, np2, all_default());
}

template <class TriangleMesh,
          class NamedParameters1>
bool
corefine_and_compute_difference(      TriangleMesh& tm1,
                                      TriangleMesh& tm2,
                                      TriangleMesh& tm_out,
                                const NamedParameters1& np1)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_difference(tm1, tm2, tm_out,
                                         np1, all_default(), all_default());
}

template <class TriangleMesh>
bool
corefine_and_compute_difference(TriangleMesh& tm1,
                                TriangleMesh& tm2,
                                TriangleMesh& tm_out)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_difference(tm1, tm2, tm_out,
                                         all_default(), all_default(), all_default());
}

///// corefine /////
template <class TriangleMesh, class NamedParameters1>
void
corefine(      TriangleMesh& tm1,
               TriangleMesh& tm2,
         const NamedParameters1& np1)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  corefine(tm1, tm2, np1, all_default());
}

template <class TriangleMesh>
void
corefine(           TriangleMesh& tm1,
                    TriangleMesh& tm2)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  corefine(tm1, tm2, all_default(), all_default());
}

#ifndef CGAL_NO_DEPRECATED_CODE
 template <class TriangleMesh,
           class NamedParameters1,
           class NamedParameters2>
 void
 corefine(      TriangleMesh& tm1,
                TriangleMesh& tm2,
          const NamedParameters1& np1,
          const NamedParameters2& np2,
          const bool throw_on_self_intersection)
{
  corefine(tm1, tm2, np1.throw_on_self_intersection(throw_on_self_intersection), np2);
}

template <class TriangleMesh, class NamedParameters1>
void
corefine(      TriangleMesh& tm1,
               TriangleMesh& tm2,
         const NamedParameters1& np1,
         const bool throw_on_self_intersection)
{
  namespace params = CGAL::Polygon_mesh_processing::parameters;
  corefine(tm1, tm2,
           np1.throw_on_self_intersection(throw_on_self_intersection),
           params::all_default());
}

template <class TriangleMesh>
void
corefine(           TriangleMesh& tm1,
                    TriangleMesh& tm2,
         const bool throw_on_self_intersection)
{
  namespace params = CGAL::Polygon_mesh_processing::parameters;
  corefine(tm1, tm2,
           params::throw_on_self_intersection(throw_on_self_intersection),
           params::all_default());
}
#endif

///// autorefine /////
namespace experimental {
template <class TriangleMesh>
void
autorefine(TriangleMesh& tm)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  autorefine(tm, all_default());
}

///// autorefine_and_remove_self_intersections /////
template <class TriangleMesh>
bool
autorefine_and_remove_self_intersections(TriangleMesh& tm)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return autorefine_and_remove_self_intersections(tm, all_default());
}

} // end of namespace experimental

} }  // end of namespace CGAL::Polygon_mesh_processing

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_COREFINEMENT_H
