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

#include <CGAL/license/Straight_skeleton_3.h>

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

namespace {

template<typename TriangleMeshIn, typename PolygonMeshOut, typename Geom_traits>
struct Output_processor
{
  using face_descriptor_in = typename boost::graph_traits<TriangleMeshIn>::face_descriptor;
  using face_descriptor_out = typename boost::graph_traits<PolygonMeshOut>::face_descriptor;

  using Polyhedron = typename CGAL::Straight_skeletons_3::internal::HDS::Polyhedron<Geom_traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using FaceGraphIO = typename CGAL::Straight_skeletons_3::IO::FaceGraphIO<Geom_traits>;

  using FT = typename Geom_traits::FT;

  // Helper to process a single result with given named parameters
  template <typename NamedParameters>
  static bool process_one(std::size_t result_index,
                          const std::vector<PolyhedronSPtr>& results_p,
                          const std::vector<FT>& save_times,
                          const TriangleMeshIn& tmesh,
                          std::vector<PolygonMeshOut>& results,
                          bool outwards,
                          const NamedParameters& np_out)
  {
    using parameters::is_default_parameter;
    using parameters::get_parameter;
    using parameters::choose_parameter;

    const FT& save_time = save_times[result_index];

    auto do_triangulate = choose_parameter(get_parameter(np_out, internal_np::do_not_triangulate_faces), false);

    std::unordered_map<face_descriptor_out, std::size_t> of2fi;
    auto of2fi_pm = boost::make_assoc_property_map(of2fi);

    CGAL_SS3_TRACE("Post processing of result @ time: " << save_time)

    PolyhedronSPtr result_p = results_p[result_index];
    CGAL_assertion(result_p && result_p->is_consistent());

    CGAL_SS3_TRACE("At time: " << save_time << ", Polyhedron with " << result_p->vertices().size()
                      << " vertices and " << result_p->facets().size() << " faces");

    PolygonMeshOut result_t;
    bool success = FaceGraphIO::save(result_p, result_t,
                                      CGAL::parameters::do_not_triangulate_faces(do_triangulate)
                                                       .face_to_face_map(of2fi_pm));
    if (!success) {
      CGAL_SS3_TRACE_V(1, "Error: failed to convert back to Surface_mesh");
      return false;
    }

    std::unordered_map<face_descriptor_out, face_descriptor_in> def_f2f;
    auto def_f2f_pm = boost::make_assoc_property_map(def_f2f);
    auto f2f = choose_parameter(get_parameter(np_out, internal_np::face_to_face_map), def_f2f_pm);

    std::cout << "f2f is provided? " << !(is_default_parameter<NamedParameters, internal_np::face_to_face_map_t>::value) << std::endl;

    for (face_descriptor_out fo : faces(result_t)) {
      std::size_t index = get(of2fi_pm, fo);
      std::cout << "index of output face " << fo << " is " << index << std::endl;
      CGAL_assertion(index < num_faces(tmesh));
      face_descriptor_in fi = *(std::next(faces(tmesh).begin(), index));
      put(f2f, fo, fi);
    }

    CGAL_postcondition(result_t.is_valid());
    CGAL_postcondition(is_valid_face_graph(result_t));
    CGAL_postcondition(CGAL::is_closed(result_t));
    CGAL_postcondition(CGAL::is_triangle_mesh(result_t));
    CGAL_postcondition(!CGAL::Polygon_mesh_processing::has_degenerate_faces(result_t));
    CGAL_postcondition(!CGAL::Polygon_mesh_processing::does_self_intersect(result_t));

    CGAL_SS3_TRACE("At time: " << save_time << ", Surface_mesh with " << num_vertices(result_t)
                      << " vertices and " << num_faces(result_t) << " faces");

    if (outwards) {
      CGAL_SS3_TRACE("Reversing face orientations...");
      CGAL::Polygon_mesh_processing::reverse_face_orientations(result_t);
    }

    results.push_back(result_t);

    return true;
  }

  // Base case: no more named parameters, use defaults for remaining results
  static bool process(std::size_t result_index,
                      const std::vector<PolyhedronSPtr>& results_p,
                      const std::vector<FT>& save_times,
                      const TriangleMeshIn& tmesh,
                      std::vector<PolygonMeshOut>& results,
                      bool outwards)
  {
    std::cout << "process(" << result_index << ", default NP)" << std::endl;

    if (result_index >= results_p.size())
      return true;

    // the variadic was shorter than the number of results, so we use default parameters for the remaining ones
    if (!process_one(result_index, results_p, save_times, tmesh, results, outwards, parameters::default_values()))
      return false;

    return process(result_index + 1, results_p, save_times, tmesh, results, outwards);
  }

  // Recursive case: process one result with the first named parameter, then recurse with remaining
  template<typename FirstNP, typename... RestNPs>
  static bool process(std::size_t result_index,
                      const std::vector<PolyhedronSPtr>& results_p,
                      const std::vector<FT>& save_times,
                      const TriangleMeshIn& tmesh,
                      std::vector<PolygonMeshOut>& results,
                      bool outwards,
                      const FirstNP& first_np,
                      const RestNPs&... rest_nps)
  {
    std::cout << "process(" << result_index << ")" << std::endl;
    std::cout << "f2f is provided? " << !(parameters::is_default_parameter<FirstNP, internal_np::face_to_face_map_t>::value) << std::endl;

    if (result_index >= results_p.size())
      return true;  // All results processed successfully

    if (!process_one(result_index, results_p, save_times, tmesh, results, outwards, first_np))
      return false;

    // the variadic was shorter than the number of results, so we use default parameters for the remaining ones
    return process(result_index + 1, results_p, save_times, tmesh, results, outwards, rest_nps...);
  }
};

} // namespace

/*!
 * \ingroup PkgStraightSkeleton3OffsettingFunctions
 *
 * This function generates face offsets of an input triangle mesh at user-specified target times.
 *
 * The offset polyhedra are constructed by progressively translating the faces of the input mesh
 * along their normal directions at a speed proportional to their weight and treating events
 * (face collapses, edge collapses, split events, etc.) as they occur. This process is
 * equivalent to the construction of a straight skeleton, and the offset polyhedra correspond
 * to the intermediate states of the input mesh during the straight skeleton construction.
 *
 * Positive time values correspond to outward offsetting, while negative values or an empty times
 * range correspond to inward offsetting. Face weights must always be positive since they represent
 * absolute speeds.
 *
 * \warning An epsilon geometric perturbation is always applied to the input mesh as to avoid
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
 *                     as key type and `FT` as value type}
 *      \cgalParamDefault{A constant property map with uniform weight `1` for all faces.}
 *      \cgalParamExtra{Precondition: all face weights must be positive.}
 *    \cgalParamNEnd
 * \cond SKIP_IN_MANUAL
 *    \cgalParamNBegin{config_file_path}
 *      \cgalParamDescription{the path to a configuration file to the algorithm. See the documentation
 *                            of the class `CGAL::Straight_skeletons_3::Configuration` for details.}
 *      \cgalParamType{`std::string`}
 *      \cgalParamDefault{no configuration file, in which case default values are used for all
 *                        configuration options (see the documentation of the class
 *                        `CGAL::Straight_skeletons_3::Configuration` for details).}
 *    \cgalParamNEnd
 * \endcond
 *  \cgalNamedParamsEnd
 *
 * \param nps_out a variadic sequence of optional \ref bgl_namedparameters "Named Parameters" to specify output options for each of the constructed offset polyhedra.
 *  \cgalNamedParamsBegin
 *   \cgalParamNBegin{do_not_triangulate_faces}
 *     \cgalParamDescription{a Boolean used to specify whether the faces of the offset meshes
 *                           should be triangulated or not.}
 *     \cgalParamDefault{`false` (i.e., faces are triangulated)}
 *     \cgalParamExtra{Note that sometimes faces must be triangulated as to be representable in
 *                     a halfedge data structure, for example faces with holes.}
 *    \cgalParamNEnd
 *    \cgalParamNBegin{face_to_face_map}
 *     \cgalParamDescription{a property map filled by this function, associating to an output face
 *                           its corresponding input face.}
 *     \cgalParamType{a class model of `WritablePropertyMap` with `std::pair<boost::graph_traits<PolygonMeshOut>::%face_descriptor, const PolygonMeshOut&>`
 *                    as key type and `boost::graph_traits<TriangleMeshIn>::%face_descriptor` as value type}
 *     \cgalParamDefault{unused}
 *   \cgalParamNEnd
 *  \cgalNamedParamsEnd
 *
 * \return `true` if offset meshes were successfully constructed; `false` otherwise.
 *
 * \pre Target times must be either all positive, or all negative.
 * \pre The input mesh is a closed, outward-oriented, triangle mesh without self-intersections.
 */
template <typename TriangleMeshIn, typename PolygonMeshOut,
          typename FT,
          typename NamedParametersIn = parameters::Default_named_parameters,
          typename... NamedParametersOut>
bool create_straight_skeleton_and_offset_polyhedra_3(const TriangleMeshIn& tmesh,
                                                     const std::vector<FT>& save_times, // intentional copy
                                                     std::vector<PolygonMeshOut>& results,
                                                     const NamedParametersIn& np_in = parameters::default_values(),
                                                     const NamedParametersOut&... nps_out)
{
  namespace SS3i = CGAL::Straight_skeletons_3::internal;
  namespace SS3io = CGAL::Straight_skeletons_3::IO;
  namespace PMP = CGAL::Polygon_mesh_processing;

  using parameters::get_parameter;
  using parameters::get_parameter_reference;
  using parameters::choose_parameter;

  using Geom_traits = typename GetGeomTraits<TriangleMeshIn, NamedParametersIn>::type;

  using Polyhedron = SS3i::HDS::Polyhedron<Geom_traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using FaceGraphIO = SS3io::FaceGraphIO<Geom_traits>;

  using face_descriptor_in = typename boost::graph_traits<TriangleMeshIn>::face_descriptor;
  using face_descriptor_out = typename boost::graph_traits<PolygonMeshOut>::face_descriptor;

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
    return false;
  }

  if (results_p.size() != save_times.size()) {
    CGAL_SS3_TRACE_V(1, "Error: expected " << save_times.size() << " polyhedra but only got " << results_p.size());
    return false;
  }

  // Start recursive processing with all named parameters
  using OutputProcessor = Output_processor<TriangleMeshIn, PolygonMeshOut, Geom_traits>;
  if (!OutputProcessor::process(0, results_p, save_times, tmesh, results, outwards, nps_out...)) {
    return false;
  }

  return true;
}

} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_CREATE_OFFSET_POLYHEDRA_H */
