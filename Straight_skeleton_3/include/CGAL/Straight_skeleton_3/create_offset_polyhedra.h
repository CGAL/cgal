// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_CREATE_OFFSET_POLYHEDRA_H
#define CGAL_STRAIGHT_SKELETON_CREATE_OFFSET_POLYHEDRA_H

#include <CGAL/Straight_skeleton_3/Configuration.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/IO/Face_graph_IO.h>
#include <CGAL/Straight_skeleton_3/create_straight_skeleton_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <filesystem>
#include <iostream>
#include <memory>
#include <vector>

namespace CGAL {

/*!
 * \ingroup PkgStraightSkeleton3OffsettingFunctions
 *
 * This function generates face offsets of an input triangle mesh at user-specified times.
 *
 * The offset polyhedra are constructed by progressively translating the faces of the input mesh
 * along their normal directions at a speed proportional to their weight and treating events
 * (face collapses, edge collapses, split events, etc.) as they occur. This process is
 * equivalent to the construction of a straight skeleton, and the offset polyhedra correspond
 * to the intermediate states of the input mesh during the straight skeleton construction.
 *
 * Positive time values correspond to outward offsetting, while negative values or an empty times range
 * correspond to inward offsetting.
 *
 * \warning An epsilon perturbation is always applied to the input mesh's geometry as to avoid
 * degenerate configurations, See \ref Straight_skeleton_3Limitations for more information.
 *
 * \tparam TriangleMeshIn must be a model of `FaceListGraph`, `HalfedgeListGraph`
 * \tparam PolygonMeshOut must be a model of `MutableFaceGraph`, `FaceListGraph`, `HalfedgeListGraph`
 * \tparam NamedParametersIn a sequence of \ref bgl_namedparameters "Named Parameters"
 * \tparam NamedParametersOut a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param tmesh the input triangle mesh whose faces are to be offset
 * \param save_times the times at which the offset polyhedra are to be constructed
 * \param results the output vector of polyhedra that will contain the constructed offset polyhedra
 * \param np_in an optional sequence of \ref bgl_namedparameters "Named Parameters"
 *              among the ones listed below
 *  \cgalNamedParamsBegin
 *    \cgalParamNBegin{vertex_point_map}
 *      \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
 *      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMeshIn>::%vertex_descriptor`
 *                     as key type and `GeomTraits::Point_3` as value type, `GeomTraits` being the type of the parameter `geom_traits`}
 *      \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
 *      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                      must be available in `TriangleMeshIn`.}
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
 *      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMeshIn>::%face_descriptor`
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
 * \param np_out an optional sequence of \ref bgl_namedparameters "Named Parameters"
 *               among the ones listed below
 *  \cgalNamedParamsBegin
 *    \cgalParamNBegin{face_weight_map}
 *     \cgalParamDescription{a property map filled by this function, associating to each output face
 *                           the weight (speed) of its corresponding face in the input.}
 *     \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<PolygonMeshOut>::%face_descriptor`
 *                    as key type and `double` as value type}
 *     \cgalParamDefault{unused}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{do_not_triangulate_faces}
 *     \cgalParamDescription{a Boolean used to specify whether the offset meshes' faces
 *                           should be triangulated or not.}
 *     \cgalParamDefault{`false` (i.e., faces are triangulated)}
 *     \cgalParamExtra{Note that sometimes faces must be triangulated as to be representable in
 *                     an halfedge data structure, for example faces with holes.}
 *    \cgalParamNEnd
 *  \cgalNamedParamsEnd
 *
 * \return `true` if offset meshes were successfully constructed; `false` otherwise.
 *
 * \pre Save offsets should be either all positive, or all negative.
 * \pre The input mesh is a closed, outward-oriented, triangle mesh without self-intersections.
 */
template <typename TriangleMeshIn, typename PolygonMeshOut,
          typename FT,
          typename NamedParametersIn = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
bool create_straight_skeleton_and_offset_polyhedra_3(const TriangleMeshIn& tmesh,
                                                     const std::vector<FT>& save_times, // intentional copy
                                                     std::vector<PolygonMeshOut>& results,
                                                     const NamedParametersIn& np_in = parameters::default_values(),
                                                     const NamedParametersOut& np_out = parameters::default_values())
{
  namespace SS3i = CGAL::Straight_skeletons_3::internal;
  namespace SS3io = CGAL::Straight_skeletons_3::IO;
  namespace PMP = CGAL::Polygon_mesh_processing;

  using Geom_traits = typename GetGeomTraits<TriangleMeshIn, NamedParametersIn>::type;

  using Polyhedron = SS3i::HDS::Polyhedron<Geom_traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using FaceGraphIO = SS3io::FaceGraphIO<Geom_traits>;

  const bool outwards = (!save_times.empty() && CGAL::is_positive(save_times.front()));

  // visitor to collect the results at save times
  struct Mesh_collector_visitor
    : public SS3i::algorithm::Default_mesh_offset_visitor<Geom_traits>
  {
    Mesh_collector_visitor(std::vector<PolyhedronSPtr>& results) : results_(results) { }

    void on_save_event(PolyhedronSPtr polyhedron, FT) override {
      PolyhedronSPtr other_polyhedron = polyhedron->clone();
      results_.push_back(other_polyhedron);
    }

  private:
    std::vector<PolyhedronSPtr>& results_;
  };

  std::vector<PolyhedronSPtr> results_p;
  Mesh_collector_visitor visitor(results_p);

  // main call
  auto ss_ptr = Straight_skeletons_3::internal::construct_skeleton(tmesh, save_times, np_in.visitor(visitor));
  if (!ss_ptr) {
    CGAL_SS3_TRACE_V(1, "Error: run() returned 'false'");
    return false;
  }

  if (results_p.size() != save_times.size()) {
    CGAL_SS3_TRACE_V(1, "Error: expected " << save_times.size() << " polyhedra but only got " << results_p.size());
    return false;
  }

  for (std::size_t i=0; i<results_p.size(); ++i) {
    const FT& save_time = save_times[i];

    CGAL_SS3_TRACE("Post processing of result @ " << save_time)

    PolyhedronSPtr result_p = results_p[i];
    CGAL_assertion(result_p && result_p->is_consistent());

    CGAL_SS3_TRACE("At time: " << save_time << ", Polyhedron with " << result_p->vertices().size()
                     << " vertices and " << result_p->facets().size() << " faces");

    // Convert back to Surface_mesh structure to save the results.
    // This could be avoided if a polygonal output was prefered, but then we
    // need some specific conversion code.
    PolygonMeshOut result_t;
    bool success = FaceGraphIO::save(result_p, result_t, np_out);
    if (!success) {
      CGAL_SS3_TRACE_V(1, "Error: failed to convert back to Surface_mesh");
      return false;
    }

    CGAL_postcondition(result_t.is_valid());
    CGAL_postcondition(is_valid_face_graph(result_t));
    CGAL_postcondition(CGAL::is_closed(result_t));
    CGAL_postcondition(CGAL::is_triangle_mesh(result_t));
    CGAL_postcondition(!PMP::has_degenerate_faces(result_t));
    CGAL_postcondition(!PMP::does_self_intersect(result_t));

    CGAL_SS3_TRACE("At time: " << save_time << ", Surface_mesh with " << num_vertices(result_t)
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

} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_CREATE_OFFSET_POLYHEDRA_H */
