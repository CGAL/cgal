// Copyright (c) 2019 GeometryFactory (France).
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
// Author(s)     : Jane Tournois

#ifndef TETRAHEDRAL_REMESHING_H
#define TETRAHEDRAL_REMESHING_H

#include <CGAL/Triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Tetrahedral_remeshing/Sizing_field.h>
#include <CGAL/Tetrahedral_remeshing/Uniform_sizing_field.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_adaptive_remeshing_impl.h>
#include <CGAL/Tetrahedral_remeshing/internal/compute_c3t3_statistics.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#ifdef CGAL_DUMP_REMESHING_STEPS
#include <sstream>
#endif

namespace CGAL
{
  ///////////////////////////////////////////////////
  ///////////////// TRIANGULATION_3 /////////////////
  ///////////////////////////////////////////////////
  /*!
  * \ingroup PkgTetrahedralRemeshingRef
  * remeshes a tetrahedral mesh.
  *
  * This function takes as input a 3-dimensional triangulation
  * and performs a sequence of atomic operations
  * in order to generate as output a high quality mesh with a prescribed density.
  * These atomic operations are performed as follows :
  *   - edge splits, until all edges satisfy a prescribed length criterion,
  *   - edge collapses, until all edges satisfy a prescribed length criterion,
  *   - edge flips, to locally improve dihedral angles, until they can't be improved by flipping,
  *   - global smoothing by vertex relocations,
  *   - re-projection of boundary vertices to the initial surface.
  *
  * This remeshing function can deal with multi-domains, multi-material boundaries and features.
  * It preserves the geometry of
  * subdomains throughout the remeshing process.
  *
  * Subdomains are defined by indices that
  * are stored in the cells of the input triangulation, following the `MeshCellBase_3`
  * concept.
  * The surfacic interfaces between subdomains are formed by facets which two incident cells
  * have different subdomain indices.
  * The edges where three or more subdomains meet form feature polylines,
  * and are considered as constrained edges.
  *
  *
  * @tparam Traits is the geometric traits, model of `RemeshingTriangulationTraits_3`
  * @tparam TDS is the triangulation data structure for `Triangulation_3`,
  *             model of ` TriangulationDataStructure_3`,
  *             with cell base model of `MeshCellBase_3`
  *             and vertex base model of `MeshVertexBase_3`.
  * @tparam SLDS is an optional parameter for `Triangulation_3`, that
  *             specifies the type of the spatial lock data structure.
  * @tparam NamedParameters a sequence of \ref Remeshing_namedparameters "Named Parameters"
  *
  * @param tr the triangulation to the remeshed, of type `Triangulation_3<Traits, TDS, SLDS>`.
  *           `Remeshing_triangulation` is a helper class that satisfies all the requirements
  *           of its template parameters.
  * @param target_edge_length the uniform target edge length. This parameter provides a
  *          mesh density target for the remeshing algorithm.
  * @param np optional sequence of \ref Remeshing_namedparameters "Named Parameters"
  *          among the ones listed below
  * \cgalNamedParamsBegin
  *  \cgalParamBegin{number_of_iterations} the number of iterations for the full
  *     sequence of atomic operations
  *     performed (listed in the above description)
  *  \cgalParamEnd
  *  \cgalParamBegin{remesh_boundaries} If `false`, none of the volume boundaries can be modified.
  *     Otherwise, the topology is preserved, but atomic operations can be performed on the
  *     surfaces, and along feature polylines, such that boundaries are remeshed.
  *  \cgalParamEnd
  *  \cgalParamBegin{edge_is_constrained_map} a property map containing the
  *    constrained - or - not status of each edge of `tr`. A constrained edge can be split
  *    or collapsed, but not flipped.
  *  \cgalParamEnd
  *  \cgalParamBegin{facet_is_constrained_map} a property map containing the
  *    constrained - or - not status of each facet of `tr`. A constrained facet can be split
  *    or collapsed, but not flipped.
  *  \cgalParamEnd
  *  \cgalParamBegin{cell_is_selected_map} a property map containing the
  *    selected - or - not status for each cell of `tr` for remeshing.
  *    Only selected cells are modified (and possibly their neighbors if surfaces are
  *    modified) by remeshing.
  *    By default, all cells with a non-zero `Subdomain_index` are selected.
  *  \cgalParamEnd
  * \cgalNamedParamsEnd

  * @todo implement 1D smoothing for constrained edges
  * @todo implement sizing field instead of uniform target edge length
  */
  template<typename Traits, typename TDS, typename SLDS,
           typename NamedParameters>
  void tetrahedral_adaptive_remeshing(
    CGAL::Triangulation_3<Traits, TDS, SLDS>& tr,
    const double& target_edge_length,
    const NamedParameters& np)
  {
    typedef CGAL::Triangulation_3<Traits, TDS, SLDS> Triangulation;
    tetrahedral_adaptive_remeshing(
      tr,
      [target_edge_length](const Triangulation::Point& p)
                          {return target_edge_length;},
      np);
  }
  
  template<typename Traits, typename TDS, typename SLDS,
           typename NamedParameters>
  void tetrahedral_adaptive_remeshing(
    CGAL::Triangulation_3<Traits, TDS, SLDS>& tr,
    const float& target_edge_length,
    const NamedParameters& np)
  {
    typedef CGAL::Triangulation_3<Traits, TDS, SLDS> Triangulation;
    tetrahedral_adaptive_remeshing(
      tr,
      [target_edge_length](const Triangulation::Point& p)
                          {return target_edge_length; },
      np);
  }

  template<typename Traits, typename TDS, typename SLDS,
           typename SizingFunction,
           typename NamedParameters>
  void tetrahedral_adaptive_remeshing(
    CGAL::Triangulation_3<Traits, TDS, SLDS>& tr,
    const SizingFunction& sizing,
    const NamedParameters& np)
  {
    CGAL_assertion(tr.is_valid(true));

    typedef CGAL::Triangulation_3<Traits, TDS, SLDS> Tr;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    bool remesh_surfaces = choose_parameter(get_parameter(np, internal_np::remesh_boundaries),
                                        true);
    bool protect = !remesh_surfaces;
    // bool adaptive = choose_parameter(get_parameter(np, internal_np::adaptive_size),
    //                              false);
    std::size_t max_it = choose_parameter(get_parameter(np, internal_np::number_of_iterations),
                                      1);

    typedef typename internal_np::Lookup_named_param_def <
      internal_np::cell_selector_t,
      NamedParameters,
      Tetrahedral_remeshing::internal::All_cells_selected<Tr>//default
    > ::type SelectionFunctor;
    SelectionFunctor cell_select
      = choose_parameter(get_parameter(np, internal_np::cell_selector),
                     Tetrahedral_remeshing::internal::All_cells_selected<Tr>());

    typedef std::pair<typename Tr::Vertex_handle, typename Tr::Vertex_handle> Edge_vv;
    typedef Tetrahedral_remeshing::internal::No_constraint_pmap<Edge_vv> No_edge;
    typedef typename internal_np::Lookup_named_param_def <
      internal_np::edge_is_constrained_t,
      NamedParameters,
      No_edge//default
    > ::type ECMap;
    ECMap ecmap = choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
                               No_edge());

    typedef typename Tr::Facet Facet;
    typedef Tetrahedral_remeshing::internal::No_constraint_pmap<Facet> No_facet;
    typedef typename internal_np::Lookup_named_param_def <
      internal_np::facet_is_constrained_t,
      NamedParameters,
      No_facet//default
    > ::type FCMap;
    FCMap fcmap = choose_parameter(get_parameter(np, internal_np::facet_is_constrained),
                               No_facet());

    typedef typename internal_np::Lookup_named_param_def <
      internal_np::remeshing_visitor_t,
      NamedParameters,
      Tetrahedral_remeshing::internal::Default_remeshing_visitor
    > ::type Visitor;
    Visitor visitor
      = choose_parameter(get_parameter(np, internal_np::remeshing_visitor),
                     Tetrahedral_remeshing::internal::Default_remeshing_visitor());

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Tetrahedral remeshing ("
      << "nb_iter = " << max_it
      << ", protect = " << std::boolalpha << protect
      << ")" << std::endl;

    std::cout << "Init tetrahedral remeshing...";
    std::cout.flush();
#endif

    typedef Tetrahedral_remeshing::internal::Adaptive_remesher<
      Tr, SizingFunction, ECMap, FCMap, SelectionFunctor, Visitor> Remesher;
    Remesher remesher(tr, sizing, protect
                    , ecmap, fcmap
                    , cell_select
                    , visitor);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "done." << std::endl;
    Tetrahedral_remeshing::internal::compute_statistics(
      remesher.tr(), cell_select, "statistics_begin.txt");
#endif

    // perform remeshing
    std::size_t nb_extra_iterations = 3;
    remesher.remesh(max_it, nb_extra_iterations);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    const double angle_bound = 5.0;
    Tetrahedral_remeshing::debug::dump_cells_with_small_dihedral_angle(tr,
      angle_bound, cell_select, "bad_cells.mesh");
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    Tetrahedral_remeshing::internal::compute_statistics(tr,
      cell_select, "statistics_end.txt");
#endif
  }

  template<typename Traits, typename TDS, typename SLDS>
  void tetrahedral_adaptive_remeshing(
    CGAL::Triangulation_3<Traits, TDS, SLDS>& tr,
    const double& target_edge_length)
  {
    tetrahedral_adaptive_remeshing(tr, target_edge_length,
      CGAL::parameters::all_default());
  }

  ///////////////////////////////////////////////////
  /////// MESH_COMPLEX_3_IN_TRIANGULATION_3 /////////
  ///////////////////////////////////////////////////

  template<typename Tr, 
           typename CornerIndex, typename CurveIndex,
           typename NamedParameters>
  void tetrahedral_adaptive_remeshing(
      CGAL::Mesh_complex_3_in_triangulation_3<Tr, CornerIndex, CurveIndex>& c3t3,
      const double& target_edge_length,
      const NamedParameters& np)
  {
    tetrahedral_adaptive_remeshing(
      c3t3,
      [target_edge_length](const typename Tr::Point& p)
                          {return target_edge_length; },
      np);
  }

  template<typename Tr,
           typename CornerIndex, typename CurveIndex,
           typename NamedParameters>
  void tetrahedral_adaptive_remeshing(
      CGAL::Mesh_complex_3_in_triangulation_3<Tr, CornerIndex, CurveIndex>& c3t3,
      const float& target_edge_length,
      const NamedParameters& np)
  {
    tetrahedral_adaptive_remeshing(
      c3t3,
      [target_edge_length](const typename Tr::Point& p)
                          {return target_edge_length; },
      np);
  }

  template<typename Tr,
           typename CornerIndex, typename CurveIndex,
           typename NamedParameters>
  void tetrahedral_adaptive_remeshing(
      CGAL::Mesh_complex_3_in_triangulation_3<Tr, CornerIndex, CurveIndex>& c3t3,
      const double& target_edge_length)
  {
    return tetrahedral_adaptive_remeshing(c3t3, target_edge_length,
                                          CGAL::parameters::all_default());
  }

  template<typename Tr,
           typename CornerIndex, typename CurveIndex,
           typename SizingFunction,
           typename NamedParameters>
  void tetrahedral_adaptive_remeshing(
      CGAL::Mesh_complex_3_in_triangulation_3<Tr, CornerIndex, CurveIndex>& c3t3,
      const SizingFunction& sizing,
      const NamedParameters& np)
  {
    CGAL_assertion(c3t3.triangulation().tds().is_valid(true));

    using parameters::get_parameter;
    using parameters::choose_parameter;

    bool remesh_surfaces = choose_parameter(get_parameter(np, internal_np::remesh_boundaries),
                                        true);
    bool protect = !remesh_surfaces;
    std::size_t max_it = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);

    typedef typename boost::lookup_named_param_def <
      internal_np::cell_selector_t,
      NamedParameters,
      Tetrahedral_remeshing::internal::All_cells_selected<Tr>//default
    > ::type SelectionFunctor;
    SelectionFunctor cell_select
      = choose_parameter(get_parameter(np, internal_np::cell_selector),
                     Tetrahedral_remeshing::internal::All_cells_selected<Tr>());

    typedef std::pair<typename Tr::Vertex_handle, typename Tr::Vertex_handle> Edge_vv;
    typedef Tetrahedral_remeshing::internal::No_constraint_pmap<Edge_vv> No_edge;
    typedef typename boost::lookup_named_param_def <
      internal_np::edge_is_constrained_t,
      NamedParameters,
      No_edge//default
    > ::type ECMap;
    ECMap ecmap = choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
                               No_edge());

    typedef typename Tr::Facet Facet;
    typedef Tetrahedral_remeshing::internal::No_constraint_pmap<Facet> No_facet;
    typedef typename boost::lookup_named_param_def <
      internal_np::facet_is_constrained_t,
      NamedParameters,
      No_facet//default
    > ::type FCMap;
    FCMap fcmap = choose_parameter(get_parameter(np, internal_np::facet_is_constrained),
                               No_facet());

    typedef typename boost::lookup_named_param_def <
      internal_np::remeshing_visitor_t,
      NamedParameters,
      Tetrahedral_remeshing::internal::Default_remeshing_visitor
    > ::type Visitor;
    Visitor visitor
      = choose_parameter(get_parameter(np, internal_np::remeshing_visitor),
                     Tetrahedral_remeshing::internal::Default_remeshing_visitor());

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Tetrahedral remeshing ("
      << "nb_iter = " << max_it
      << "protect = " << std::boolalpha << protect << ", "
      << ")" << std::endl;

    std::cout << "Init tetrahedral remeshing...";
    std::cout.flush();
#endif

    typedef Tetrahedral_remeshing::internal::Adaptive_remesher<
      Tr, SizingFunction, ECMap, FCMap, SelectionFunctor,
      Visitor,
      CornerIndex, CurveIndex
    > Remesher;
    Remesher remesher(c3t3, sizing, protect
                    , ecmap, fcmap
                    , cell_select
                    , visitor);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "done." << std::endl;
    Tetrahedral_remeshing::internal::compute_statistics(
      remesher.tr(),
      cell_select, "statistics_begin.txt");
#endif

    // perform remeshing
    std::size_t nb_extra_iterations = 3;
    remesher.remesh(max_it, nb_extra_iterations);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    const double angle_bound = 5.0;
    Tetrahedral_remeshing::debug::dump_cells_with_small_dihedral_angle(
      c3t3.triangulation(),
      angle_bound, cell_select, "bad_cells.mesh");
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    Tetrahedral_remeshing::internal::compute_statistics(
      c3t3.triangulation(),
      cell_select, "statistics_end.txt");
#endif
  }


}//end namespace CGAL

#endif //TETRAHEDRAL_REMESHING_H
