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

#include <CGAL/Straight_skeleton_3/Configuration.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/IO/Face_graph_IO.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Straight_skeleton_builder_3.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_perturbation.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/Named_function_parameters.h>

#include <memory>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {

template <typename TriangleMesh,
          typename FT,
          typename NamedParameters>
std::shared_ptr<CGAL::Straight_skeleton_3<typename GetGeomTraits<TriangleMesh, NamedParameters>::type> >
construct_skeleton(const TriangleMesh& tmesh,
                   std::vector<FT> save_times, // intentional copy
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

  using Polyhedron = SS3i::HDS::Polyhedron<Geom_traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Perturbation = SS3i::algorithm::Polyhedron_perturbation<Geom_traits>;
  using Straight_skeleton_builder_3 = SS3i::algorithm::Straight_skeleton_builder_3<Geom_traits>;
  using FaceGraphIO = SS3io::FaceGraphIO<Geom_traits>;

  // Get config file path from named parameters if provided, else use default
  SS3::ConfigurationSPtr config = SS3::Configuration::get_instance();
  config->load_default_values();

  if constexpr (!parameters::is_default_parameter<NamedParameters, internal_np::config_file_path_t>::value) {
    std::string str_conf_file = parameters::get_parameter(np, internal_np::config_file_path);
    CGAL_SS3_IO_TRACE("Loading configuration from file: " << str_conf_file);
    if (!config->load(str_conf_file)) {
      CGAL_SS3_IO_TRACE("Error: Failed to load configuration file: " << str_conf_file);
      return { };
    }
  } else {
    CGAL_SS3_IO_TRACE("Using default configuration values (no config file path provided)");
  }

  const std::filesystem::path save_path = choose_parameter(get_parameter(np, internal_np::io_path),
                                                           std::filesystem::current_path());

  CGAL_precondition(!CGAL::is_empty(tmesh));
  CGAL_precondition(CGAL::is_closed(tmesh));
  CGAL_precondition(!PMP::does_self_intersect(tmesh));
  CGAL_precondition(PMP::is_outward_oriented(tmesh));

  const bool outwards = (!save_times.empty() && CGAL::is_positive(save_times.front()));

  // check that all times are of the same sign
  for (const FT& time : save_times) {
    if (CGAL::is_zero(time)) {
      CGAL_SS3_TRACE_V(1, "Error: save time should be non-zero");
      return { };
    } else if (CGAL::is_positive(time) != outwards) {
      CGAL_SS3_TRACE_V(1, "Error: save times must be all positive or all negative.");
      return { };
    }
  }

  // The underlying implementation always shrinks the mesh, so for outwards offsetting,
  // we reverse the face orientations (whether the mesh was outward oriented or not).
  if (outwards) {
    for (FT& t : save_times) {
      t = -t;
    }
  }

  // Convert the suface mesh into the SLS3-specific data structure that allows faces with multiple
  // borders and disconnected facet connected components
  PolyhedronSPtr p = FaceGraphIO::convert(tmesh, np.outward_offsetting(outwards));
  CGAL_SS3_DEBUG_SPTR(p);
  CGAL_SS3_TRACE("Post conversion: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  Perturbation::apply_rand_perturbation(p);
  CGAL_SS3_TRACE("Post perturbation: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  using Default_visitor = SS3i::algorithm::Default_mesh_offset_visitor<Geom_traits>;
  using Visitor = typename internal_np::Lookup_named_param_def<internal_np::visitor_t,
                                                               NamedParameters,
                                                               Default_visitor>::reference;

  Default_visitor default_visitor;
  Visitor visitor = choose_parameter(get_parameter_reference(np, internal_np::visitor), default_visitor);

  auto skel_builder = Straight_skeleton_builder_3::create(p, save_times, save_path);
  skel_builder->set_visitor(&visitor);
  skel_builder->set_outward(outwards);

  bool success = skel_builder->run();
  if (!success) {
    CGAL_SS3_TRACE_V(1, "Error: run() returned 'false'");
    return { };
  }

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("final_skeleton.obj", skel_builder->get_skeleton());
#endif

  return skel_builder->get_skeleton();
}

} // namespace internal
} // namespace Straight_skeletons_3

/*!
 * \ingroup PkgStraightSkeleton3SkeletonFunctions
 *
 * constructs the interior or exterior straight skeleton of a 3D polyhedron.
 *
 * Positive time values signify outward construction, while negative values or an empty times range
 * correspond to inward construction.
 *
 * Face weights may be passed to represent different front speeds. These values must always be positive
 * since they represent absolute speeds.
 *
 * The upper bound `max_time` defines how far the straight skeleton is constructed,
 * i.e. nodes at a time greater than `max_time` are not created. This may result
 * in the presence of unbounded arcs and sheets in the straight skeleton.
 *
 * \warning An epsilon geometric perturbation is always applied to the input mesh as to avoid
 * degenerate configurations, See \ref Straight_skeleton_3Limitations for more information.
 *
 * \tparam TriangleMesh must be a model of `FaceListGraph`, `HalfedgeListGraph`
 * \tparam FT must be a model of `FieldNumberType` compatible with the type of `geom_traits`
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
                        where `GeomTraits` being the type of the parameter `geom_traits`}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{geom_traits}
 *      \cgalParamDescription{an instance of a geometric traits class}
 *      \cgalParamType{a class model of `Kernel`}
 *      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type
 *                      and must provide exact predicates and exact constructions.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{face_weight_map}
 *      \cgalParamDescription{a property map associating to each face the weight (speed) of the face.}
 *      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *                     as key type and `FT` as value type}
 *      \cgalParamDefault{A constant property map with uniform weight `1` for all faces.}
 *      \cgalParamExtra{Precondition: all face weights must be positive.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{config_file_path}
 *      \cgalParamDescription{the path to a configuration file to the algorithm. See the documentation
 *                            of the class `CGAL::Straight_skeletons_3::Configuration` for details.}
 *      \cgalParamType{`std::string`}
 *      \cgalParamDefault{A set of default values, see the documentation of the configuration class.}
 *    \cgalParamNEnd
 *  \cgalNamedParamsEnd
 *
 * \pre The value `max_time` must be either positive or negative (non-zero).
 * \pre The input `tmesh` must be a non-empty, closed, and self-intersection-free triangle mesh.
 */
template <typename TriangleMesh,
          typename FT,
          typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
std::shared_ptr<CGAL::Straight_skeleton_3<geom_traits> >
#else
std::shared_ptr<CGAL::Straight_skeleton_3<typename GetGeomTraits<TriangleMesh, NamedParameters>::type> >
#endif
create_straight_skeleton_3(const TriangleMesh& tmesh,
                           const FT& max_time,
                           const NamedParameters& np = parameters::default_values())
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
 */
template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
std::shared_ptr<CGAL::Straight_skeleton_3<geom_traits> >
#else
std::shared_ptr<CGAL::Straight_skeleton_3<typename GetGeomTraits<TriangleMesh, NamedParameters>::type> >
#endif
create_straight_skeleton_3(const TriangleMesh& tmesh,
                           const NamedParameters& np = parameters::default_values())
{
  using FT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type::FT;
  return Straight_skeletons_3::internal::construct_skeleton(tmesh, std::vector<FT>{ }, np);
}

} // namespace CGAL

#endif // CGAL_CREATE_STRAIGHT_SKELETON_3_H
