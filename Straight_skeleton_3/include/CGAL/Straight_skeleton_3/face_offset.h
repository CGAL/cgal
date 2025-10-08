// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_FACE_OFFSET_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_FACE_OFFSET_H

#include <CGAL/Straight_skeleton_3/Configuration.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/IO/Face_graph_IO.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Straight_skeleton_builder.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_perturbation.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_self_intersection.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <filesystem>
#include <iostream>
#include <memory>
#include <vector>

namespace CGAL {
namespace Straight_skeletons_3 {

/*!
 * \ingroup PkgStraightSkeleton3OffsettingFunctions
 *
 * Given a range of values, this function generates face offsets of an input triangle mesh.
 *
 * The offset meshes are constructed by translating the faces of the input mesh by a distance
 * that depends on the specified offset and their respective weights (speeds). During
 * offsetting, elements of the mesh interact with each other (merging, splitting, etc.).
 *
 * Positive offset values correspond to outward offsets, while negative values correspond to inward offsets.
 *
 * \warning An epsilon perturbation is always applied to the input mesh's geometry as to avoid so-called
 * degenerate positions, meaning a configuration where any three supporting planes of the facets
 * of the input mesh do not intersect in a single point.
 *
 * \tparam TriangleMesh must be a model of `FaceListGraph`, `HalfedgeListGraph`
 * \tparam PolygonMeshOut must be a model of `MutableFaceGraph`, `FaceListGraph`, `HalfedgeListGraph`
 * \tparam NamedParametersIn must be a model of `NamedParameters`
 * \tparam NamedParametersOut must be a model of `NamedParameters`
 *
 * \param tmesh the input triangle mesh whose faces are to be offset
 * \param save_times the times at which the offset meshes are to be constructed
 * \param results the output vector of meshes that will contain the constructed offset meshes
 * \param np_in the input named parameters
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
 *    \cgalParamNBegin{face_weight}
 *      \cgalParamDescription{a property map associating to each face the weight (speed) of the face.}
 *      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *                     as key type and `double` as value type}
 *      \cgalParamDefault{A constant property map with uniform weight 1.0 for all faces.}
 *      \cgalParamExtra{Precondition: all face weights must be positive.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{config_file_path}
 *      \cgalParamDescription{the path to a configuration file to the algorithm. See the documentation
 *                            of the `Configuration` class for details.}
 *      \cgalParamType{`std::string`}
 *      \cgalParamDefault{The path to a default configuration file which must be found in the working directory.}
 *    \cgalParamNEnd
 *  \cgalNamedParamsEnd
 *
 * \param np_out the output named parameters
 *  \cgalNamedParamsBegin
 *    \cgalParamNBegin{face_weight}
 *      \cgalParamDescription{a property map filled by this function, associating to each face
 *                            the weight (speed) of the corresponding face in the input.}
 *      \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
 *                     as key type and `double` as value type}
 *      \cgalParamDefault{unused}
 *    \cgalParamNEnd
 *  \cgalNamedParamsEnd
 *
 * \pre save offsets should be either all positive, or all negative.
 * \pre the input mesh is a closed, outward-oriented, triangle mesh without self-intersections.
 *
 * \return `true` if offset meshes were successfully constructed; `false` otherwise.
 */
template <typename TriangleMeshIn, typename PolygonMeshOut,
          typename FT,
          typename NamedParametersIn = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
bool face_offset(TriangleMeshIn& tmesh,
                 std::vector<FT> save_times, // intentional copy
                 std::vector<PolygonMeshOut>& results,
                 const NamedParametersIn& np_in = parameters::default_values(),
                 const NamedParametersOut& np_out = parameters::default_values())
{
  namespace SS3i = CGAL::Straight_skeletons_3::internal;
  namespace PMP = CGAL::Polygon_mesh_processing;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  using Geom_traits = typename GetGeomTraits<TriangleMeshIn, NamedParametersIn>::type;

  using Polyhedron = SS3i::HDS::Polyhedron<Geom_traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Transformation = SS3i::algorithm::PolyhedronTransformation<Geom_traits>;
  using Perturbation = SS3i::algorithm::PolyhedronPerturbation<Geom_traits>;
  using SelfIntersection = SS3i::algorithm::SelfIntersection<Geom_traits>;
  using SimpleStraightSkel = SS3i::algorithm::SimpleStraightSkel<Geom_traits>;
  using FaceGraphIO = IO::FaceGraphIO<Geom_traits>;

  // Get config file path from named parameters if provided, else use default (working directory)
  ConfigurationSPtr config = Configuration::getInstance();
  std::string str_conf_file = parameters::choose_parameter(parameters::get_parameter(np_in, internal_np::config_file_path),
                                                           config->findDefaultFilename());

  std::cout << "Seek config file @ " << str_conf_file << std::endl;
  if (!config->load(str_conf_file)) {
    std::cerr << "Error: Config file '" << str_conf_file << "' not found." << std::endl;
    return false;
  }

  const std::filesystem::path save_path = choose_parameter(get_parameter(np_in, internal_np::io_path),
                                                           std::filesystem::current_path());

  CGAL_precondition(!CGAL::is_empty(tmesh));
  CGAL_precondition(CGAL::is_closed(tmesh));
  CGAL_precondition(!CGAL::is_triangle_mesh(tmesh) || PMP::is_outward_oriented(tmesh));
  CGAL_precondition(!CGAL::is_triangle_mesh(tmesh) || !PMP::does_self_intersect(tmesh));

  const bool outwards = (!save_times.empty() && CGAL::is_positive(save_times.front()));

  // check that all offsets are of the same sign (ignoring zeros)
  for (const FT& time : save_times) {
    if (CGAL::is_zero(time)) {
      std::cerr << "Error: time should be non-zero" << std::endl;
      return false;
    } else if (CGAL::is_positive(time) != outwards) {
      std::cerr << "Error: offsets must all be positive or all negative." << std::endl;
      return false;
    }
  }

  // The underlying implementation always shrinks the mesh, so for outwards offsets,
  // we reverse the face orientations, whether the mesh was outward oriented or not.
  if (outwards) {
    CGAL_SS3_TRACE("Reversing face orientations...");
    PMP::reverse_face_orientations(tmesh);

    // also switch to negative offsets
    for (FT& t : save_times) {
      t = -t;
    }
  }

  // Convert the suface mesh into the SLS3-specific data structure that allows faces with multiple
  // borders and disconnected facet connected components
  PolyhedronSPtr p = FaceGraphIO::convert(tmesh, np_in);
  CGAL_SS3_DEBUG_SPTR(p);
  Transformation::normalizeFacetPlanes(p);

  CGAL_SS3_TRACE("Post conversion: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  // Perturbation to ensure generic configuration
  bool safe_mode = true;
  if (config->isLoaded()) {
    if ((config->contains("Preprocessing", "check_degenerate_configuration") &&
         !config->getBool("Preprocessing", "check_degenerate_configuration"))) {
      safe_mode = false;
    }
  }

  PolyhedronSPtr p_mem = p->clone();

  // We always need to ensure that points are exactly on the planes of their incident facets
  Perturbation::randTiltPlanesv3(p);

  if (safe_mode) {
    for (;;) {
      if (Perturbation::doAll2PlanesIntersect(p) &&
          Perturbation::doAll3PlanesIntersect(p) &&
          !SelfIntersection::hasSelfIntersectingSurface(p)) {
        CGAL_SS3_TRACE("Found a good perturbation");
        break;
      }

      p = p_mem->clone();
      Perturbation::randTiltPlanesv3(p);
    }
  }

  CGAL_SS3_TRACE("Post perturbation: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  // visitor to collect the results at save offsets
  struct Mesh_collector_visitor
    : public internal::algorithm::Default_mesh_offset_visitor<Geom_traits>
  {
    Mesh_collector_visitor(std::vector<PolyhedronSPtr>& results) : results_(results) { }

    void on_save_event(PolyhedronSPtr polyhedron, FT) override {
      PolyhedronSPtr other_polyhedron = polyhedron->clone();
      results_.push_back(other_polyhedron);
    }

  private:
    std::vector<PolyhedronSPtr>& results_;
  };

  auto skel_builder = SimpleStraightSkel::create(p, save_times, save_path);

  std::vector<PolyhedronSPtr> results_p;
  Mesh_collector_visitor visitor(results_p);
  skel_builder->setVisitor(&visitor);

  // main call
  bool success = skel_builder->run();
  if (!success) {
    return false;
  }

  if (results_p.size() != save_times.size()) {
    return false;
  }

  if (outwards) {
    for (FT& t : save_times) {
      t = -t;
    }
  }

  for (std::size_t i=0; i<results_p.size(); ++i) {
    const FT& save_time = save_times[i];

    CGAL_SS3_TRACE("Post processing of result @ " << save_time)

    PolyhedronSPtr result_p = results_p[i];
    CGAL_assertion(result_p && result_p->isConsistent());

    CGAL_SS3_TRACE("At offset: " << save_time << ", Polyhedron with " << result_p->vertices().size()
                     << " vertices and " << result_p->facets().size() << " faces");

    // Convert back to Surface_mesh structure to save the results.
    // This could be avoided if a polygonal output was prefered, but then we
    // need some specific conversion code.
    PolygonMeshOut result_t;
    bool success = IO::FaceGraphIO<Geom_traits>::save(result_p, result_t, np_out);
    if (!success) {
      std::cerr << "Error: failed to convert back to Surface_mesh" << std::endl;
      return false;
    }

    CGAL_postcondition(result_t.is_valid());
    CGAL_postcondition(is_valid_face_graph(result_t));
    CGAL_postcondition(CGAL::is_closed(result_t));
    CGAL_postcondition(CGAL::is_triangle_mesh(result_t));
    CGAL_postcondition(!PMP::has_degenerate_faces(result_t));
    CGAL_postcondition(!PMP::does_self_intersect(result_t));

    CGAL_SS3_TRACE("At offset: " << save_time << ", Surface_mesh with " << num_vertices(result_t)
                     << " vertices and " << num_faces(result_t) << " faces");

    if (outwards) {
      CGAL_SS3_TRACE("Reversing face orientations...");
      PMP::reverse_face_orientations(result_t);
    }

    results.push_back(result_t);
  }

  CGAL_SS3_TRACE("Offset mesh(es) constructed");
  return true;
}

} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_FACE_OFFSET_H */
