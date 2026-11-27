// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_CREATE_STRAIGHT_SKELETON_3_H
#define CGAL_CREATE_STRAIGHT_SKELETON_3_H

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/internal/Configuration.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/IO/Face_graph_IO.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Straight_skeleton_builder_3.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/Named_function_parameters.h>

#include <filesystem>
#include <memory>
#include <type_traits>
#include <vector>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {

template <typename TriangleMesh,
          typename IFT,
          typename NamedParameters>
std::shared_ptr<CGAL::Straight_skeleton_3<typename GetGeomTraits<TriangleMesh, NamedParameters>::type> >
construct_skeleton(const TriangleMesh& tmesh,
                   const std::vector<IFT>& input_save_times,
                   const NamedParameters& np)
{
  namespace SS3 = CGAL::Straight_skeletons_3;
  namespace SS3i = CGAL::Straight_skeletons_3::internal;
  namespace SS3io = CGAL::Straight_skeletons_3::IO;
  namespace PMP = CGAL::Polygon_mesh_processing;

  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;
  using parameters::get_parameter_reference;

  using Geom_traits = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using FT = typename Geom_traits::FT;

  static_assert(std::is_same_v<typename CGAL::Algebraic_structure_traits<FT>::Is_exact, CGAL::Tag_true>,
                "The construction of 3D straight skeletons requires a kernel with exact constructions.");

  using Polyhedron = SS3i::HDS::Polyhedron<Geom_traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Perturbation = SS3i::algorithm::Polyhedron_perturbation<Geom_traits>;
  using Straight_skeleton_builder_3 = SS3i::algorithm::Straight_skeleton_builder_3<Geom_traits>;
  using FaceGraphIO = SS3io::FaceGraphIO<Geom_traits>;

  // Get config file path from named parameters if provided, else use default
  SS3::ConfigurationSPtr config = SS3::Configuration::get_instance();
  config->load_default_values();

  if constexpr (!parameters::is_default_parameter<NamedParameters, internal_np::config_file_path_t>::value) {
    std::filesystem::path config_path = parameters::get_parameter(np, internal_np::config_file_path);
    if (std::filesystem::exists(config_path) && std::filesystem::is_regular_file(config_path)) {
      CGAL_SS3_IO_TRACE_V(1, "Loading configuration from file: " << config_path);
      if (!config->load(config_path)) {
        CGAL_SS3_IO_TRACE_V(1, "Error: Failed to load configuration file: " << config_path);
        return { };
      }
    } else {
      CGAL_SS3_IO_TRACE_V(1, "Warning: provided config path leads nowhere: " << config_path);
    }
  } else {
    CGAL_SS3_IO_TRACE_V(1, "Using default configuration values (no config file path provided)");
  }

  const std::filesystem::path save_path = choose_parameter(get_parameter(np, internal_np::io_path),
                                                           std::filesystem::current_path());

  CGAL_precondition_code(auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                                      get_const_property_map(vertex_point, tmesh));)
  CGAL_precondition_code(auto pmp_np = parameters::vertex_point_map(vpm).geom_traits(
                                         choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits)));)
  CGAL_precondition(!CGAL::is_empty(tmesh));
  CGAL_precondition(CGAL::is_closed(tmesh));
  CGAL_warning(CGAL::is_triangle_mesh(tmesh));
  CGAL_precondition(!CGAL::is_triangle_mesh(tmesh) || !PMP::does_self_intersect(tmesh, pmp_np));
  CGAL_precondition(!CGAL::is_triangle_mesh(tmesh) || PMP::is_outward_oriented(tmesh, pmp_np));

  const std::vector<FT> save_times(input_save_times.begin(), input_save_times.end());

  for (const FT& time : save_times) {
    if (!CGAL::is_positive(time)) {
      CGAL_SS3_TRACE_V(1, "Error: save times must be positive");
      return { };
    }
  }

  // The sign of the weights is the direction of the propagation: faces move along their normal
  // (outward) if their weight is positive, and against it (inward) if their weight is negative.
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  auto weight_pmap = choose_parameter(get_parameter(np, internal_np::face_weight),
                                      CGAL::Constant_property_map<face_descriptor, FT>(-1));

  CGAL::Sign weight_sign = CGAL::ZERO;
  for (const face_descriptor f : faces(tmesh)) {
    const FT weight = get(weight_pmap, f);
    const CGAL::Sign s = CGAL::sign(weight);
    if (s == CGAL::ZERO) {
      CGAL_SS3_TRACE_V(1, "Error: face weights must be non-zero");
      return { };
    } else if (weight_sign != CGAL::ZERO && s != weight_sign) {
      CGAL_SS3_TRACE_V(1, "Error: face weights must be all positive or all negative");
      return { };
    }
    weight_sign = s;
  }

  // an exterior straight skeleton is unbounded, and is only constructed up to a given time
  if (weight_sign == CGAL::POSITIVE && save_times.empty()) {
    CGAL_SS3_TRACE_V(1, "Error: a maximum time is required when face weights are positive");
    return { };
  }

  // Convert the suface mesh into the SPS3-specific data structure that allows faces with multiple
  // borders and disconnected facet connected components
  PolyhedronSPtr p = FaceGraphIO::convert(tmesh, np);
  CGAL_SS3_DEBUG_SPTR(p);
  CGAL_SS3_TRACE_V(2, "Post conversion: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  Perturbation::apply_rand_perturbation(p);
  CGAL_SS3_TRACE_V(2, "Post perturbation: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  using Default_visitor = SS3i::algorithm::Default_mesh_offset_visitor<Geom_traits>;
  using Visitor = typename internal_np::Lookup_named_param_def<internal_np::visitor_t,
                                                               NamedParameters,
                                                               Default_visitor>::reference;

  Default_visitor default_visitor;
  Visitor visitor = choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor);

  auto skel_builder = Straight_skeleton_builder_3::create(p, save_times, save_path);
  skel_builder->set_visitor(&visitor);

  bool success = skel_builder->run();
  if (!success) {
    CGAL_SS3_TRACE_V(1, "Error: run() returned 'false'");
    return { };
  }

  return skel_builder->get_skeleton();
}

} // namespace internal
} // namespace Straight_skeletons_3

/*!
 * \ingroup PkgStraightSkeleton3SkeletonFunctions
 *
 * constructs the interior or exterior straight skeleton of a 3D polyhedron.
 *
 * The weight of a face is its signed speed along its oriented normal
 * (see Section \ref Straight_skeleton_3WeightedStraightSkeletonDef): faces move inward if weights
 * are negative, which constructs the interior skeleton, and outward if weights are positive,
 * which constructs the exterior skeleton. In both cases, times are positive: time `0` corresponds
 * to the input polyhedron.
 *
 * The bound `max_time` defines how far the straight skeleton is constructed,
 * i.e. nodes at a time greater than `max_time` are not computed. This may result
 * in the presence of unbounded spokes and sheets in the straight skeleton.
 *
 * \warning An epsilon geometric perturbation is always applied to the input mesh so as to avoid
 * degenerate configurations; see \ref Straight_skeleton_3GeneralPosition for more information.
 *
 * \tparam TriangleMesh must be a model of `FaceListGraph`, `HalfedgeListGraph`
 * \tparam FT must be a number type convertible to `geom_traits::FT`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param tmesh the input triangle mesh whose straight skeleton is to be constructed
 * \param max_time the maximum time until which the straight skeleton is to be constructed
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters"
 *           among the ones listed below
 *  \cgalNamedParamsBegin
 *    \cgalParamNBegin{vertex_point_map}
 *      \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
 *      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                     as key type and a model of `Kernel::Point_3` as value type}
 *      \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
 *      \cgalParamExtra{The value type of the vertex point property map must be compatible with `GeomTraits::Point_3`,
 *                      where `GeomTraits` being the type of the parameter `geom_traits`}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{geom_traits}
 *      \cgalParamDescription{an instance of a geometric traits class}
 *      \cgalParamType{a class model of `Kernel`}
 *      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type
 *                      and must provide exact predicates and exact constructions.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{face_weight_map}
 *      \cgalParamDescription{a property map associating to each face its weight (speed), that is the signed speed
 *                            of the face along its outward normal.}
 *      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *                     as key type and `FT` as value type}
 *      \cgalParamDefault{A constant property map with uniform weight `-1` for all faces.}
 *      \cgalParamExtra{Precondition: face weights must be non-zero, and must all have the same sign.}
 *      \cgalParamExtra{Adjacent coplanar faces must have the same weight.}
 *    \cgalParamNEnd
 * \cond SKIP_IN_MANUAL
 *    \cgalParamNBegin{config_file_path}
 *      \cgalParamDescription{the path to a configuration file to the algorithm. See the documentation
 *                            of the class `CGAL::Straight_skeletons_3::Configuration` for details.}
 *      \cgalParamType{`std::filesystem::path`}
 *      \cgalParamDefault{no configuration file, in which case default values are used for all
 *                        configuration options (see the documentation of the class
 *                        `CGAL::Straight_skeletons_3::Configuration` for details).}
 *    \cgalParamNEnd
 * \endcond
 *  \cgalNamedParamsEnd
 *
 * \returns a pointer to the straight skeleton, or a null pointer if the construction failed.
 *
 * \pre The value `max_time` must be positive.
 * \pre The input `tmesh` must be a non-empty, closed, outward-oriented, and self-intersection-free triangle mesh.
 */
template <typename TriangleMesh,
          typename FT,
          typename CGAL_NP_TEMPLATE_PARAMETERS
#ifndef DOXYGEN_RUNNING
          , typename = std::enable_if_t<!is_named_function_parameter<FT> >
#endif
          >
#ifdef DOXYGEN_RUNNING
std::shared_ptr<CGAL::Straight_skeleton_3<geom_traits> >
#else
std::shared_ptr<CGAL::Straight_skeleton_3<typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type> >
#endif
create_straight_skeleton_3(const TriangleMesh& tmesh,
                           const FT& max_time,
                           const CGAL_NP_CLASS& np = parameters::default_values())
{
  // this indirection is because we add a visitor when we are constructing offset polyhedra
  return Straight_skeletons_3::internal::construct_skeleton(tmesh, std::vector<FT>{ max_time }, np);
}

/*!
 * \ingroup PkgStraightSkeleton3SkeletonFunctions
 *
 * constructs the complete interior straight skeleton of a 3D polyhedron.
 *
 * See the other overload for a comprehensive documentation of the parameters.
 *
 * \returns a pointer to the interior straight skeleton, or a null pointer if the construction failed.
 *
 * \pre All face weights must be negative.
 * \pre The input `tmesh` must be a non-empty, closed, outward-oriented, and self-intersection-free triangle mesh.
 */
template <typename TriangleMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
std::shared_ptr<CGAL::Straight_skeleton_3<geom_traits> >
#else
std::shared_ptr<CGAL::Straight_skeleton_3<typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type> >
#endif
create_straight_skeleton_3(const TriangleMesh& tmesh,
                           const CGAL_NP_CLASS& np = parameters::default_values())
{
  using FT = typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT;
  return Straight_skeletons_3::internal::construct_skeleton(tmesh, std::vector<FT>{ }, np);
}

} // namespace CGAL

#endif // CGAL_CREATE_STRAIGHT_SKELETON_3_H
