// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_FACET_OFFSET_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_FACET_OFFSET_H

#include <CGAL/Straight_skeleton_3/Configuration.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/IO/Surface_meshIO.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Straight_skeleton_builder.h>

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
* Given a range of values, this function generates face-offsets of an input mesh.
*
* The mesh is constructed by translating the faces of the input mesh by a distance
* that depends on the specified offset and their respective weights (speeds); during
* offsetting, elements of the mesh interact with each other (merging, splitting, etc.).
*
* Positive offset values correspond to outward offsets, while negative values correspond to inward offsets.
*
* \warning A tiny perturbation is always applied to the input mesh's geometry as to avoid so-called
* degenerate positions, meaning a configuration where any three supporting planes of the facets
* of the input mesh do not intersect in a single point.
*
* \tparam TriangleMesh must be a model of `MutableFaceGraph`, `FaceListGraph`, `HalfedgeListGraph`
* \tparam NamedParametersIn must be a model of `NamedParameters`
* \tparam NamedParametersOut must be a model of `NamedParameters`
*
* \param tmesh the input mesh to offset
* \param save_offsets the offsets to apply to the mesh
* \param results the output vector of meshes that will contain the offset meshes
* \param np_in the input named parameters
* \param np_out the output named parameters
*
* \pre offsets value should be non-zero and all of the same sign.
* \return `true` if offset meshes were successfully constructed; `false` otherwise.
*/
template <typename TriangleMesh,
          typename FT,
          typename NamedParametersIn = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
bool face_offset(TriangleMesh& tmesh,
                 std::vector<FT> save_offsets, // intentional copy
                 std::vector<TriangleMesh>& results,
                 const NamedParametersIn& np_in = parameters::default_values(),
                 const NamedParametersOut& np_out = parameters::default_values())
{
  namespace SS3i = CGAL::Straight_skeletons_3::internal;
  namespace PMP = CGAL::Polygon_mesh_processing;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  using Geom_traits = typename GetGeomTraits<TriangleMesh, NamedParametersIn>::type;

  using Polyhedron = SS3i::HDS::Polyhedron<Geom_traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using PolyhedronTransformation = SS3i::algorithm::PolyhedronTransformation<Geom_traits>;
  using SimpleStraightSkel = SS3i::algorithm::SimpleStraightSkel<Geom_traits>;
  using Surface_meshIO = IO::Surface_meshIO<Geom_traits>;

  // @tmp move to somewhere sensible
  // SLS3 config file
  ConfigurationSPtr config = Configuration::getInstance();
  std::string str_conf_file = config->findDefaultFilename();
  if (!config->load(str_conf_file)) {
    std::cerr << "Error: Config file '" << str_conf_file << "' not found." << std::endl;
    return false;
  }

  const std::filesystem::path save_path = choose_parameter(get_parameter(np_in, internal_np::io_path),
                                                           std::filesystem::current_path());

  CGAL_precondition(!CGAL::is_empty(tmesh));
  CGAL_precondition(CGAL::is_closed(tmesh));
  CGAL_precondition(CGAL::is_triangle_mesh(tmesh));
  CGAL_precondition(!PMP::does_self_intersect(tmesh));

  const bool outwards = (!save_offsets.empty() && CGAL::is_positive(save_offsets.front()));

  // The underlying implementation always shrinks the mesh, so for outwards offsets,
  // we reverse the face orientations, whether the mesh was outward oriented or not.
  if (outwards) {
    CGAL_SS3_TRACE("Reversing face orientations...");
    PMP::reverse_face_orientations(tmesh);

    // also switch to negative offsets
    for (FT& offset : save_offsets) {
      offset = -offset;
    }
  }

  // Convert the suface mesh into the SLS3-specific data structure that allows faces with multiple
  // borders and disconnected facet connected components
  PolyhedronSPtr p = Surface_meshIO::convert(tmesh, np_in);

  CGAL_SS3_TRACE("Post merge: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  PolyhedronTransformation::truncatePrecision(p);

  // @todo In safe mode, this should be:
  // if (bad) -> reduce amplitude of perturbation and try again"
  PolyhedronTransformation::normalizeFacetPlanes(p);
  PolyhedronTransformation::randTiltPlanesv3(p);
  CGAL_assertion(p->isConsistent());

  CGAL_SS3_TRACE("Post randomization: " << p->vertices().size() << " NV " << p->facets().size() << " NF");

  // visitor to collect the results
  struct Mesh_collector_visitor
    : public internal::algorithm::Default_mesh_offset_visitor<Geom_traits>
  {
    Mesh_collector_visitor(std::vector<PolyhedronSPtr>& results) : results_(results) { }

    void on_save_offset_event(PolyhedronSPtr polyhedron, FT) override {
      PolyhedronSPtr other_polyhedron = polyhedron->clone();
      results_.push_back(other_polyhedron);
    }

  private:
    std::vector<PolyhedronSPtr>& results_;
  };

  // run the skeleton code
  auto skel_builder = SimpleStraightSkel::create(p, save_offsets, save_path);

  std::vector<PolyhedronSPtr> results_p;
  Mesh_collector_visitor visitor(results_p);
  skel_builder->setVisitor(&visitor);

  bool success = skel_builder->run();
  if (!success) {
    return false;
  }

  if (results_p.size() != save_offsets.size()) {
    return false;
  }

  if (outwards) {
    for (FT& offset : save_offsets) {
      offset = -offset;
    }
  }

  for (std::size_t i=0; i<results_p.size(); ++i) {
    const FT& save_offset = save_offsets[i];

    CGAL_SS3_TRACE("Post processing of result @ " << save_offset)

    PolyhedronSPtr result_p = results_p[i];
    CGAL_assertion(result_p && result_p->isConsistent());

    CGAL_SS3_TRACE("At offset: " << save_offset << ", Polyhedron with " << result_p->vertices().size()
                     << " vertices and " << result_p->facets().size() << " faces");

    // Convert back to Surface_mesh structure to save the results.
    // This could be avoided if a polygonal output was prefered, but then we
    // need some specific conversion code.
    TriangleMesh result_t;
    bool success = IO::Surface_meshIO<Geom_traits>::save(result_p, result_t, np_out);
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

    CGAL_SS3_TRACE("At offset: " << save_offset << ", Surface_mesh with " << num_vertices(result_t)
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

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_FACET_OFFSET_H */
