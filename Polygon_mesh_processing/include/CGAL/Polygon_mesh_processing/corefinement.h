// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_COREFINEMENT_H
#define CGAL_POLYGON_MESH_PROCESSING_COREFINEMENT_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Visitor.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Face_graph_output_builder.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/Output_builder_for_autorefinement.h>
#include <CGAL/iterator.h>

namespace CGAL {

#if !defined(CGAL_NO_DEPRECATED_CODE) && !defined(DOXYGEN_RUNNING)
namespace Corefinement {
using Polygon_mesh_processing::Corefinement::Self_intersection_exception;
}
#endif

namespace Polygon_mesh_processing {

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


#define CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP(I) \
  typedef typename internal_np::Lookup_named_param_def < \
    internal_np::edge_is_constrained_t, \
    NamedParametersOut##I, \
    Corefinement::No_mark<TriangleMesh> \
  > ::type Ecm_out_##I; \
    Ecm_out_##I ecm_out_##I = \
      parameters::choose_parameter<Ecm_out_##I>(parameters::get_parameter(std::get<I>(nps_out), internal_np::edge_is_constrained));

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
  * @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters"
  * @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters"
  * @tparam NamedParametersOut0 a sequence of \ref bgl_namedparameters "Named Parameters" for computing the union of the volumes bounded by `tm1` and `tm2`
  * @tparam NamedParametersOut1 a sequence of \ref bgl_namedparameters "Named Parameters" for computing the intersection of the volumes bounded by `tm1` and `tm2`
  * @tparam NamedParametersOut2 a sequence of \ref bgl_namedparameters "Named Parameters" for computing the difference of the volumes bounded by `tm1` and `tm2`
  * @tparam NamedParametersOut3 a sequence of \ref bgl_namedparameters "Named Parameters" for computing the difference of the volumes bounded by `tm2` and `tm1`
  *
  * @param tm1 first input triangulated surface mesh
  * @param tm2 second input triangulated surface mesh
  * @param output an array of output surface meshes
  * @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  * @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm1` (`tm2`)}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm1 (tm2))`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
  *                     must be available in `TriangleMesh`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{edge_is_constrained_map}
  *     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm1` (`tm2`)}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
  *                    as key type and `bool` as value type}
  *     \cgalParamDefault{a constant property map returning `false` for any edge}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{face_index_map}
  *     \cgalParamDescription{a property map associating to each face of `tm1` (`tm2`) a unique index
  *                           between `0` and `num_faces(tm1) - 1` (`num_faces(tm2) - 1`)}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
  *                    as key type and `std::size_t` as value type}
  *     \cgalParamDefault{an automatically indexed internal map}
  *     \cgalParamExtra{If the property map is writable, the indices of the faces of `tm1` and `tm2`
  *                     will be set after the corefinement is done.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{visitor}
  *     \cgalParamDescription{a visitor used to track the creation of new faces}
  *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
  *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
  *     \cgalParamExtra{`np1` only}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{throw_on_self_intersection}
  *     \cgalParamDescription{If `true`, the set of triangles closed to the intersection of `tm1` and `tm2` will be
  *                           checked for self-intersections and `Corefinement::Self_intersection_exception`
  *                           will be thrown if at least one self-intersection is found.}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *     \cgalParamExtra{`np1` only}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @param nps_out an optional tuple of sequences of \ref bgl_namedparameters "Named Parameters" each among the ones listed below
  *        (`tm_out` being used to refer to the output surface mesh in `output` corresponding to a given named parameter sequence)
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm_out`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm_out)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
  *                     must be available in `TriangleMesh`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{edge_is_constrained_map}
  *     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm_out`.
  *                           An edge of `tm_out` is constrained if it is on the intersection of `tm1` and `tm2`,
  *                           or if the edge corresponds to a constrained edge in `tm1` or `tm2`.}
  *     \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
  *                    as key type and `bool` as value type}
  *   \cgalParamNEnd
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
std::array<bool,4>
corefine_and_compute_boolean_operations(
        TriangleMesh& tm1,
        TriangleMesh& tm2,
  const std::array< boost::optional<TriangleMesh*>,4>& output,
  const NamedParameters1& np1,
  const NamedParameters2& np2,
  const std::tuple<NamedParametersOut0,
                     NamedParametersOut1,
                     NamedParametersOut2,
                     NamedParametersOut3>& nps_out)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  const bool throw_on_self_intersection =
    choose_parameter(get_parameter(np1, internal_np::throw_on_self_intersection), false);

// Vertex point maps
  //for input meshes
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters1>::type  VPM1;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters2>::type  VPM2;

  CGAL_static_assertion((std::is_same<typename boost::property_traits<VPM1>::value_type,
                                      typename boost::property_traits<VPM2>::value_type>::value));

  VPM1 vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                               get_property_map(boost::vertex_point, tm1));

  VPM2 vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                               get_property_map(boost::vertex_point, tm2));

  typedef typename boost::property_traits<VPM1>::value_type                 Point_3;

  // for output meshes: here we have to use a trick so that if for a specific output
  // that is not requested, the default vpm does not have the same value type as the
  // input map, a dummy default vpm is used so that calls to get/put can be compiled
  // (even if not used).
  typedef std::tuple<
    Corefinement::TweakedGetVertexPointMap<Point_3, NamedParametersOut0, TriangleMesh>,
    Corefinement::TweakedGetVertexPointMap<Point_3, NamedParametersOut1, TriangleMesh>,
    Corefinement::TweakedGetVertexPointMap<Point_3, NamedParametersOut2, TriangleMesh>,
    Corefinement::TweakedGetVertexPointMap<Point_3, NamedParametersOut3, TriangleMesh>
  > VPM_out_tuple_helper;

  typedef std::tuple<
    boost::optional< typename std::tuple_element<0, VPM_out_tuple_helper>::type::type >,
    boost::optional< typename std::tuple_element<1, VPM_out_tuple_helper>::type::type >,
    boost::optional< typename std::tuple_element<2, VPM_out_tuple_helper>::type::type >,
    boost::optional< typename std::tuple_element<3, VPM_out_tuple_helper>::type::type >
  > VPM_out_tuple;

  VPM_out_tuple vpm_out_tuple(
    Corefinement::get_vpm<Point_3>(std::get<0>(nps_out), output[0],
                                   typename std::tuple_element<0, VPM_out_tuple_helper>::type::Use_default_tag()),
    Corefinement::get_vpm<Point_3>(std::get<1>(nps_out), output[1],
                                   typename std::tuple_element<1, VPM_out_tuple_helper>::type::Use_default_tag()),
    Corefinement::get_vpm<Point_3>(std::get<2>(nps_out), output[2],
                                   typename std::tuple_element<2, VPM_out_tuple_helper>::type::Use_default_tag()),
    Corefinement::get_vpm<Point_3>(std::get<3>(nps_out), output[3],
                                   typename std::tuple_element<3, VPM_out_tuple_helper>::type::Use_default_tag())
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
                        parameters::vertex_point_map(*std::get<Corefinement::UNION>(vpm_out_tuple)));
    if (output[Corefinement::INTERSECTION] != boost::none)
      if (&tm1 != *output[Corefinement::INTERSECTION])
        copy_face_graph(tm1,
                        *(*output[Corefinement::INTERSECTION]),
                        parameters::vertex_point_map(vpm1),
                        parameters::vertex_point_map(*std::get<Corefinement::INTERSECTION>(vpm_out_tuple)));

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
                        parameters::vertex_point_map(*std::get<Corefinement::UNION>(vpm_out_tuple)));
    if (output[Corefinement::INTERSECTION] != boost::none)
      clear(*(*output[Corefinement::INTERSECTION]));
    if (output[Corefinement::TM1_MINUS_TM2] != boost::none)
      clear(*(*output[Corefinement::TM1_MINUS_TM2]));
    if (output[Corefinement::TM2_MINUS_TM1] != boost::none)
      if (&tm2 != *output[Corefinement::TM2_MINUS_TM1])
        copy_face_graph(tm2,
                        *(*output[Corefinement::TM2_MINUS_TM1]),
                        parameters::vertex_point_map(vpm2),
                        parameters::vertex_point_map(*std::get<Corefinement::TM2_MINUS_TM1>(vpm_out_tuple)));
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
                          parameters::vertex_point_map(*std::get<Corefinement::UNION>(vpm_out_tuple)));
      if (output[Corefinement::INTERSECTION] != boost::none)
        clear(*(*output[Corefinement::INTERSECTION]));
      if (output[Corefinement::TM2_MINUS_TM1] != boost::none)
        clear(*(*output[Corefinement::TM2_MINUS_TM1]));
      if (output[Corefinement::TM1_MINUS_TM2] != boost::none)
        if (&tm1 != *output[Corefinement::TM1_MINUS_TM2])
          copy_face_graph(tm1,
                          *(*output[Corefinement::TM1_MINUS_TM2]),
                          parameters::vertex_point_map(vpm1),
                          parameters::vertex_point_map(*std::get<Corefinement::TM1_MINUS_TM2>(vpm_out_tuple)));
      return CGAL::make_array(true, true, true, true);
    }

// Edge is-constrained maps
  //for input meshes
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters1,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm1;

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters2,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm2;

  Ecm1 ecm1 = choose_parameter<Ecm1>(get_parameter(np1, internal_np::edge_is_constrained));
  Ecm2 ecm2 = choose_parameter<Ecm2>(get_parameter(np2, internal_np::edge_is_constrained));

  typedef Corefinement::Ecm_bind<TriangleMesh, Ecm1, Ecm2> Ecm_in;

  //for output meshes
  CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP(0)
  CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP(1)
  CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP(2)
  CGAL_COREF_SET_OUTPUT_EDGE_MARK_MAP(3)

  // In the current version all types must be the same so an array would be fine
  typedef std::tuple<Ecm_out_0, Ecm_out_1, Ecm_out_2, Ecm_out_3>
                                                            Edge_mark_map_tuple;

  // Face index point maps
  typedef typename CGAL::GetInitializedFaceIndexMap<TriangleMesh, NamedParameters1>::type FaceIndexMap1;
  typedef typename CGAL::GetInitializedFaceIndexMap<TriangleMesh, NamedParameters2>::type FaceIndexMap2;

  FaceIndexMap1 fid_map1 = get_initialized_face_index_map(tm1, np1);
  FaceIndexMap2 fid_map2 = get_initialized_face_index_map(tm2, np2);


  // User visitor
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::visitor_t,
    NamedParameters1,
    Corefinement::Default_visitor<TriangleMesh>//default
  > ::type User_visitor;
  User_visitor uv(choose_parameter<User_visitor>(get_parameter(np1, internal_np::visitor)));

  // surface intersection algorithm call
  typedef Corefinement::Face_graph_output_builder<TriangleMesh,
                                                  VPM1,
                                                  VPM2,
                                                  VPM_out_tuple,
                                                  FaceIndexMap1,
                                                  FaceIndexMap2,
                                                  Default,
                                                  Ecm_in,
                                                  Edge_mark_map_tuple,
                                                  User_visitor> Ob;

  typedef Corefinement::Surface_intersection_visitor_for_corefinement<
            TriangleMesh, VPM1, VPM2, Ob, Ecm_in, User_visitor> Algo_visitor;

  Ecm_in ecm_in(tm1,tm2,ecm1,ecm2);
  Edge_mark_map_tuple ecms_out(ecm_out_0, ecm_out_1, ecm_out_2, ecm_out_3);
  Ob ob(tm1, tm2, vpm1, vpm2, fid_map1, fid_map2, ecm_in, vpm_out_tuple, ecms_out, uv, output);

  // special case used for clipping open meshes
  if (choose_parameter(get_parameter(np1, internal_np::use_bool_op_to_clip_surface), false))
  {
    CGAL_assertion(output[Corefinement::INTERSECTION] != boost::none);
    CGAL_assertion(output[Corefinement::UNION] == boost::none);
    CGAL_assertion(output[Corefinement::TM1_MINUS_TM2] == boost::none);
    CGAL_assertion(output[Corefinement::TM2_MINUS_TM1] == boost::none);
    const bool use_compact_clipper =
      choose_parameter(get_parameter(np1, internal_np::use_compact_clipper), true);
    ob.setup_for_clipping_a_surface(use_compact_clipper);
  }

  Corefinement::Intersection_of_triangle_meshes<TriangleMesh, VPM1, VPM2, Algo_visitor >
    functor(tm1, tm2, vpm1, vpm2, Algo_visitor(uv,ob,ecm_in));
  functor(CGAL::Emptyset_iterator(), throw_on_self_intersection, true);


  return CGAL::make_array(ob.union_is_valid(),
                          ob.intersection_is_valid(),
                          ob.tm1_minus_tm2_is_valid(),
                          ob.tm2_minus_tm1_is_valid());
}

template <class TriangleMesh>
std::array<bool,4>
corefine_and_compute_boolean_operations(
        TriangleMesh& tm1,
        TriangleMesh& tm2,
  const std::array< boost::optional<TriangleMesh*>,4>& output)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_boolean_operations(tm1, tm2, output,
                                                 all_default(), all_default(),
                                                 std::make_tuple(all_default(), all_default(),
                                                                   all_default(), all_default()));
}

template <class TriangleMesh, class NamedParameters1>
std::array<bool,4>
corefine_and_compute_boolean_operations(
        TriangleMesh& tm1,
        TriangleMesh& tm2,
  const std::array< boost::optional<TriangleMesh*>,4>& output,
  const NamedParameters1& np1)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_boolean_operations(tm1, tm2, output,
                                                 np1, all_default(),
                                                 std::make_tuple(all_default(), all_default(),
                                                                   all_default(), all_default()));
}

template <class TriangleMesh, class NamedParameters1, class NamedParameters2>
std::array<bool,4>
corefine_and_compute_boolean_operations(
        TriangleMesh& tm1,
        TriangleMesh& tm2,
  const std::array< boost::optional<TriangleMesh*>,4>& output,
  const NamedParameters1& np1,
  const NamedParameters2& np2)
{
  using namespace CGAL::Polygon_mesh_processing::parameters;
  return corefine_and_compute_boolean_operations(tm1, tm2, output,
                                                 np1, np2,
                                                 std::make_tuple(all_default(), all_default(),
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
  * @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters"
  * @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters"
  * @tparam NamedParametersOut a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param tm1 first input triangulated surface mesh
  * @param tm2 second input triangulated surface mesh
  * @param tm_out output surface mesh
  * @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  * @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm1` (`tm2`)}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm1 (tm2))`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
  *                     must be available in `TriangleMesh`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{edge_is_constrained_map}
  *     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm1` (`tm2`)}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
  *                    as key type and `bool` as value type}
  *     \cgalParamDefault{a constant property map returning `false` for any edge}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{face_index_map}
  *     \cgalParamDescription{a property map associating to each face of `tm1` (`tm2`) a unique index
  *                           between `0` and `num_faces(tm1) - 1` (`num_faces(tm2) - 1`)}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
  *                    as key type and `std::size_t` as value type}
  *     \cgalParamDefault{an automatically indexed internal map}
  *     \cgalParamExtra{If the property map is writable, the indices of the faces of `tm1` and `tm2`
  *                     will be set after the corefinement is done.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{visitor}
  *     \cgalParamDescription{a visitor used to track the creation of new faces}
  *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
  *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
  *     \cgalParamExtra{`np1` only}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{throw_on_self_intersection}
  *     \cgalParamDescription{If `true` the set of triangles closed to the intersection of `tm1` and `tm2` will be
  *                           checked for self-intersections and `Corefinement::Self_intersection_exception`
  *                           will be thrown if at least one self-intersection is found.}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *     \cgalParamExtra{`np1` only}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @param np_out an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tm_out`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm_out)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
  *                     must be available in `TriangleMesh`.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{edge_is_constrained_map}
  *     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm_out`.
  *                           An edge of `tm_out` is constrained if it is on the intersection of `tm1` and `tm2`,
  *                           or if the edge corresponds to a constrained edge in `tm1` or `tm2`.}
  *     \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
  *                    as key type and `bool` as value type}
  *   \cgalParamNEnd
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
  std::array< boost::optional<TriangleMesh*>,4> output;
  output[Corefinement::UNION]=&tm_out;

  return
   corefine_and_compute_boolean_operations(tm1, tm2, output, np1, np2,
                                           std::make_tuple(np_out,
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
  std::array< boost::optional<TriangleMesh*>,4> output;
  output[Corefinement::INTERSECTION]=&tm_out;

  return
    corefine_and_compute_boolean_operations(tm1, tm2, output, np1, np2,
                                            std::make_tuple(all_default(),
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
  std::array< boost::optional<TriangleMesh*>,4> output;
  output[TM1_MINUS_TM2]=&tm_out;

  return
    corefine_and_compute_boolean_operations(tm1, tm2, output, np1, np2,
                                            std::make_tuple(all_default(),
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
 * @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters"
 * @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tm1 first input triangulated surface mesh
 * @param tm2 second input triangulated surface mesh
 * @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 * @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm1` (`tm2`)}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm1 (tm2))`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{edge_is_constrained_map}
 *     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm1` (`tm2`)}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
 *                    as key type and `bool` as value type}
 *     \cgalParamDefault{a constant property map returning `false` for any edge}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{visitor}
 *     \cgalParamDescription{a visitor used to track the creation of new faces}
 *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
 *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
 *     \cgalParamExtra{`np1` only}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{throw_on_self_intersection}
 *     \cgalParamDescription{If `true` the set of triangles closed to the intersection of `tm1` and `tm2` will be
 *                           checked for self-intersections and `Corefinement::Self_intersection_exception`
 *                           will be thrown if at least one self-intersection is found.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *     \cgalParamExtra{`np1` only}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{do_not_modify}
 *     \cgalParamDescription{if `true`, the corresponding mesh will not be updated.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *     \cgalParamExtra{If this parameter is set to `true` for both meshes nothing will be done.
 *                      If this option is set to `true` for one mesh,
 *                      the other mesh is no longer required to be without self-intersection.}
 *   \cgalParamNEnd
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
  using parameters::choose_parameter;
  using parameters::get_parameter;

  TriangleMesh* const_mesh_ptr=nullptr;
  if (choose_parameter(get_parameter(np1, internal_np::do_not_modify), false))
  {
    if (choose_parameter(get_parameter(np2, internal_np::do_not_modify), false))
      return;
    const_mesh_ptr=&tm1;
  }
  else
  {
    if (choose_parameter(get_parameter(np2, internal_np::do_not_modify), false))
      const_mesh_ptr=&tm2;
  }

  const bool throw_on_self_intersection =
    choose_parameter(get_parameter(np1, internal_np::throw_on_self_intersection), false);

// Vertex point maps
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters1>::type VPM1;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters2>::type VPM2;

  CGAL_static_assertion((std::is_same<typename boost::property_traits<VPM1>::value_type,
                                      typename boost::property_traits<VPM2>::value_type>::value));

  VPM1 vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                               get_property_map(boost::vertex_point, tm1));

  VPM2 vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                               get_property_map(boost::vertex_point, tm2));

// Edge is-constrained maps
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters1,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm1;

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters2,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm2;

  Ecm1 ecm1 = choose_parameter<Ecm1>(get_parameter(np1, internal_np::edge_is_constrained));
  Ecm2 ecm2 = choose_parameter<Ecm2>(get_parameter(np2, internal_np::edge_is_constrained));

  typedef Corefinement::Ecm_bind<TriangleMesh, Ecm1, Ecm2> Ecm;

  if (&tm1==&tm2)
  {
    Corefinement::mark_all_edges(tm1, ecm1);
    Corefinement::mark_all_edges(tm2, ecm2);
    return;
  }

  // User visitor
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::visitor_t,
    NamedParameters1,
    Corefinement::Default_visitor<TriangleMesh>//default
  > ::type User_visitor;
  User_visitor uv(choose_parameter<User_visitor>(get_parameter(np1, internal_np::visitor)));

// surface intersection algorithm call
  typedef Corefinement::No_extra_output_from_corefinement<TriangleMesh> Ob;
  typedef Corefinement::Surface_intersection_visitor_for_corefinement<
    TriangleMesh, VPM1, VPM2, Ob, Ecm, User_visitor> Algo_visitor;

  Ob ob;
  Ecm ecm(tm1,tm2,ecm1,ecm2);
  Corefinement::Intersection_of_triangle_meshes<TriangleMesh, VPM1, VPM2, Algo_visitor>
    functor(tm1, tm2, vpm1, vpm2, Algo_visitor(uv,ob,ecm,const_mesh_ptr));
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
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *  \cgalParamNEnd
 *
 *   \cgalParamNBegin{edge_is_constrained_map}
 *     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
 *                    as key type and `bool` as value type}
 *     \cgalParamDefault{a constant property map returning `false` for any edge}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{face_index_map}
 *     \cgalParamDescription{a property map associating to each face of `tm` a unique index between `0` and `num_faces(tm) - 1`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{an automatically indexed internal map}
 *     \cgalParamExtra{If the property map is writable, the indices of the faces of `tm1` and `tm2`
 *                     will be set after the corefinement is done.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{visitor}
 *     \cgalParamDescription{a visitor used to track the creation of new faces}
 *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
 *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh,
          class NamedParameters>
void
autorefine(      TriangleMesh& tm,
           const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

// Vertex point maps
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VPM;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(boost::vertex_point, tm));

// Edge is-constrained maps
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm;
  Ecm ecm = choose_parameter<Ecm>(get_parameter(np, internal_np::edge_is_constrained));

// User visitor
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::visitor_t,
    NamedParameters,
    Corefinement::Default_visitor<TriangleMesh>//default
  > ::type User_visitor;
  User_visitor uv(choose_parameter<User_visitor>(get_parameter(np, internal_np::visitor)));


// surface intersection algorithm call
  typedef Corefinement::No_extra_output_from_corefinement<TriangleMesh> Ob;
  typedef Corefinement::Surface_intersection_visitor_for_corefinement<
    TriangleMesh, VPM, VPM, Ob, Ecm, User_visitor,true> Algo_visitor;
  Ob ob;

  Corefinement::Intersection_of_triangle_meshes<TriangleMesh, VPM, VPM, Algo_visitor>
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
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{edge_is_constrained_map}
 *     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
 *                    as key type and `bool` as value type}
 *     \cgalParamDefault{a constant property map returning `false` for any edge}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{face_index_map}
 *     \cgalParamDescription{a property map associating to each face of `tm` a unique index between `0` and `num_faces(tm) - 1`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *                    as key type and `std::size_t` as value type}
 *     \cgalParamDefault{an automatically indexed internal map}
 *     \cgalParamExtra{If the property map is writable, the indices of the faces of `tm` will be set
 *                     after the autorefinement is done.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{visitor}
 *     \cgalParamDescription{a visitor used to track the creation of new faces}
 *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
 *     \cgalParamDefault{`Corefinement::Default_visitor<TriangleMesh>`}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh,
          class NamedParameters>
bool
autorefine_and_remove_self_intersections(      TriangleMesh& tm,
                                         const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

// Vertex point maps
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(boost::vertex_point, tm));

// Face index map
  typedef typename GetInitializedFaceIndexMap<TriangleMesh, NamedParameters>::type Fid_map;
  Fid_map fid_map = get_initialized_face_index_map(tm, np);


// Edge is-constrained maps
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters,
    Corefinement::No_mark<TriangleMesh>//default
  > ::type Ecm;
  Ecm ecm = choose_parameter<Ecm>(get_parameter(np, internal_np::edge_is_constrained));

// User visitor
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::visitor_t,
    NamedParameters,
    Corefinement::Default_visitor<TriangleMesh>//default
  > ::type User_visitor;
  User_visitor uv(choose_parameter<User_visitor>(get_parameter(np, internal_np::visitor)));

// surface intersection algorithm call
  typedef Corefinement::Output_builder_for_autorefinement<TriangleMesh,
                                                          VPM,
                                                          Fid_map,
                                                          Ecm,
                                                          Default > Ob;

  typedef Corefinement::Surface_intersection_visitor_for_corefinement<
    TriangleMesh, VPM, VPM, Ob, Ecm, User_visitor,true> Algo_visitor;
  Ob ob(tm, vpm, fid_map, ecm);

  Corefinement::Intersection_of_triangle_meshes<TriangleMesh, VPM, VPM, Algo_visitor>
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
