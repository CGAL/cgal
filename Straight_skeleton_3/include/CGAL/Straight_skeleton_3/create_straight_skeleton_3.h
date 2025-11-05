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

#include <CGAL/Straight_skeleton_3/Configuration.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/IO/Face_graph_IO.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Straight_skeleton_builder_3.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_perturbation.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_self_intersection.h>

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

  using parameters::get_parameter;
  using parameters::get_parameter_reference;
  using parameters::choose_parameter;

  using Geom_traits = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;

  using Polyhedron = SS3i::HDS::Polyhedron<Geom_traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Transformation = SS3i::algorithm::Polyhedron_transformation<Geom_traits>;
  using Perturbation = SS3i::algorithm::Polyhedron_perturbation<Geom_traits>;
  using Self_intersection = SS3i::algorithm::Self_intersection<Geom_traits>;
  using Straight_skeleton_builder_3 = SS3i::algorithm::Straight_skeleton_builder_3<Geom_traits>;
  using FaceGraphIO = SS3io::FaceGraphIO<Geom_traits>;

  // Get config file path from named parameters if provided, else use default (working directory)
  SS3::ConfigurationSPtr config = SS3::Configuration::get_instance();
  std::string str_conf_file = parameters::choose_parameter(parameters::get_parameter(np, internal_np::config_file_path),
                                                           config->find_default_filename());

  CGAL_SS3_TRACE_V(4, "Seek config file @ " << str_conf_file);
  if (!config->load(str_conf_file)) {
    CGAL_SS3_TRACE_V(1, "Error: Config file '" << str_conf_file << "' not found.");
    return { };
  }

  const std::filesystem::path save_path = choose_parameter(get_parameter(np, internal_np::io_path),
                                                           std::filesystem::current_path());

  CGAL_precondition(!CGAL::is_empty(tmesh));
  CGAL_precondition(CGAL::is_closed(tmesh));
  CGAL_precondition(!PMP::does_self_intersect(tmesh));
  CGAL_precondition(PMP::is_outward_oriented(tmesh));

  const bool outwards = (!save_times.empty() && CGAL::is_positive(save_times.front()));

  // check that all times are of the same sign (ignoring zeros)
  for (const FT& time : save_times) {
    if (CGAL::is_zero(time)) {
      CGAL_SS3_TRACE_V(1, "Error: time should be non-zero");
      return { };
    } else if (CGAL::is_positive(time) != outwards) {
      CGAL_SS3_TRACE_V(1, "Error: times must all be positive or all negative.");
      return { };
    }
  }

  // The underlying implementation always shrinks the mesh, so for outwards offsetting,
  // we reverse the face orientations, whether the mesh was outward oriented or not.
  if (outwards) {
    for (FT& t : save_times) {
      t = -t;
    }
  }

  // Convert the suface mesh into the SLS3-specific data structure that allows faces with multiple
  // borders and disconnected facet connected components
  PolyhedronSPtr p = FaceGraphIO::convert(tmesh, np.outward_offsetting(outwards));
  CGAL_SS3_DEBUG_SPTR(p);
  Transformation::normalize_facet_planes(p);

  CGAL_SS3_TRACE("Post conversion: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  // Perturbation to ensure generic configuration
  bool safe_mode = true;
  if (config->is_loaded()) {
    if ((config->contains("Preprocessing", "check_degenerate_configuration") &&
         !config->get_Boolean("Preprocessing", "check_degenerate_configuration"))) {
      safe_mode = false;
    }
  }

  PolyhedronSPtr p_mem = p->clone();

  // We always need to ensure that points are exactly on the planes of their incident facets
  Perturbation::apply_rand_plane_tilts_V3(p);

  if (safe_mode) {
    for (;;) {
      if (Perturbation::do_all_plane_pairs_intersect(p) &&
          Perturbation::do_all_plane_triplets_intersect(p) &&
          !Self_intersection::has_self_intersecting_surface(p)) {
        CGAL_SS3_TRACE("Found a good perturbation");
        break;
      }

      p = p_mem->clone();
      Perturbation::apply_rand_plane_tilts_V3(p);
    }
  }

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
 *                     as key type and `GeomTraits::Point_3` as value type, `GeomTraits` being the type of the parameter `geom_traits`}
 *      \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
 *      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                      must be available in `TriangleMesh`.}
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
 *                     as key type and `double` as value type}
 *      \cgalParamDefault{A constant property map with uniform weight 1.0 for all faces.}
 *      \cgalParamExtra{Precondition: all face weights must be positive.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{config_file_path}
 *      \cgalParamDescription{the path to a configuration file to the algorithm. See the documentation
 *                            of the class `CGAL::Straight_skeletons_3::Configuration` for details.}
 *      \cgalParamType{`std::string`}
 *      \cgalParamDefault{The path to a default configuration file which must be found in the working directory.}
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
