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


#include <CGAL/Tetrahedral_remeshing/Sizing_field.h>
#include <CGAL/Tetrahedral_remeshing/Uniform_sizing_field.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_adaptive_remeshing_impl.h>
#include <CGAL/Tetrahedral_remeshing/internal/compute_c3t3_statistics.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#ifdef CGAL_DUMP_REMESHING_STEPS
#include <sstream>
#endif

namespace CGAL
{
  /*!
  * \ingroup PkgTetrahedralRemeshingRef
  * remeshes a tetrahedral mesh.
  *
  * This function takes as input a 3-dimensional triangulation
  * and performs a sequence of atomic operations
  * in order to generate as output a high quality mesh with a prescribed density.
  * These atomic operations are performed as follows :
  *   - edge splits, until all edges satisfy a prescribed length criterion,
  *   - edge collapses, ntil all edges satisfy a prescribed length criterion,
  *   - edge flips, to locally improve dihedral angles, until they can't be improved by flipping,
  *   - global smoothing by vertex relocations,
  *   - re-projection of boundary vertices to the initial surface.
  *
  * This remeshing function can deal with multi-domains, multi-material boundaries and features.
  * It preserves the geometry of
  * subdomains throughout the remeshing process.
  *
  * Subdomains are defined by indices that
  * are stored in the cells of the input triangulation, following the `RemeshingCellBase_3`
  * concept.
  * The surfacic interfaces between subdomains are formed by facets which two incident cells
  * have different subdomain indices.
  * The edges where three or more subdomains meet form feature polylines,
  * and are considered as constrained edges.
  *
  *
  * @tparam Triangulation a 3-dimensional triangulation
  * deriving from `Triangulation_3`,
  * with geometric traits model of `RemeshingTriangulationTraits_3`,
  * cell base model of `RemeshingCellBase_3`
  * and vertex base model of `RemeshingVertexBase_3`.
  *
  * @tparam NamedParameters a sequence of \ref Remeshing_namedparameters "Named Parameters"
  *
  * @param tr the triangulation to the remeshed
  * @param target_edge_length the uniform target edge length. This parameter provides a
  *          mesh density target for the remeshing algorithm.
  * @param np optional sequence of \ref Remeshing_namedparameters "Named Parameters"
  *          among the ones listed below
  * \cgalNamedParamsBegin
  *  \cgalParamBegin{number_of_iterations} the number of iterations for the full
  *     sequence of atomic operations
  *     performed (listed in the above description)
  *  \cgalParamEnd
  *  \cgalParamBegin{protect_boundaries} If `true`, none of the volume boundaries can be modified.
  *     Otherwise, the geometry is preserved, but atomic operations can be performed on the
  *     surfaces, and along feature polylines.
  *  \cgalParamEnd
  *  \cgalParamBegin{edge_is_constrained_map} a property map containing the
  *    constrained - or - not status of each edge of `tr`. A constrained edge can be split
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

  //  * @tparam SizingField model of `CGAL::Sizing_field`

  //*  \cgalParamBegin{ adaptive } If `true`, size of elements adapts
  //*     ....
  //*  \cgalParamEnd
    //template<typename Triangulation,
  //         typename SizingField,
  //         typename NamedParameters>
  //void tetrahedral_adaptive_remeshing(Triangulation& tr,
  //                                    const SizingField& sizing_field,
  //                                    const NamedParameters& np)
  template<typename Triangulation, typename NamedParameters>
  void tetrahedral_adaptive_remeshing(Triangulation& tr,
                                      const double& target_edge_length,
                                      const NamedParameters& np)
  {
    tetrahedral_adaptive_remeshing(
      tr,
      [target_edge_length](const typename Triangulation::Point& p)
                          {return target_edge_length;},
      np);
  }

  template<typename Triangulation, typename NamedParameters>
  void tetrahedral_adaptive_remeshing(Triangulation& tr,
    const float& target_edge_length,
    const NamedParameters& np)
  {
    tetrahedral_adaptive_remeshing(
      tr,
      [target_edge_length](const typename Triangulation::Point& p)
                          {return target_edge_length; },
      np);
  }

  template<typename Triangulation,
           typename SizingFunction,
           typename NamedParameters>
  void tetrahedral_adaptive_remeshing(Triangulation& tr,
                                      const SizingFunction& sizing,
                                      const NamedParameters& np)
  {
    CGAL_assertion(tr.is_valid(true));

    typedef Triangulation Tr;

    using boost::choose_param;
    using boost::get_param;

    bool protect = choose_param(get_param(np, internal_np::protect_boundaries),
                                false);
    // bool adaptive = choose_param(get_param(np, internal_np::adaptive_size),
    //                              false);
    std::size_t max_it = choose_param(get_param(np, internal_np::number_of_iterations),
                                      1);

    typedef typename boost::lookup_named_param_def <
      internal_np::cell_selector_t,
      NamedParameters,
      Tetrahedral_remeshing::internal::All_cells_selected<Tr>//default
    > ::type SelectionFunctor;
    SelectionFunctor cell_select
      = choose_param(get_param(np, internal_np::cell_selector),
                     Tetrahedral_remeshing::internal::All_cells_selected<Tr>());

    typedef std::pair<typename Tr::Vertex_handle, typename Tr::Vertex_handle> Edge_vv;
    typedef Tetrahedral_remeshing::internal::No_constraint_pmap<Edge_vv> No_constraint;

    typedef typename boost::lookup_named_param_def <
      internal_np::edge_is_constrained_t,
      NamedParameters,
      No_constraint//default
    > ::type ECMap;
    ECMap ecmap = choose_param(get_param(np, internal_np::edge_is_constrained)
                             , No_constraint());

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Tetrahedral remeshing ("
      << "nb_iter = " << max_it
      << "protect = " << std::boolalpha << protect << ", "
      << ")" << std::endl;

    std::cout << "Init tetrahedral remeshing...";
    std::cout.flush();
#endif

    typedef Tetrahedral_remeshing::internal::Adaptive_remesher<
      Tr, SizingFunction, ECMap, SelectionFunctor> Remesher;
    Remesher remesher(tr, sizing, protect, ecmap
                    , cell_select
                  /*, adaptive*/);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "done." << std::endl;
    Tetrahedral_remeshing::internal::compute_statistics(
      remesher.triangulation(),
      remesher.imaginary_index(), cell_select, "statistics_begin.txt");
#endif

    remesher.preprocess();

    std::size_t it_nb = 0;
    while (it_nb++ < max_it)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "# Iteration " << it_nb << " #" << std::endl;
#endif
      if (!remesher.resolution_reached())
      {
        remesher.split();
        remesher.collapse();
      }
      remesher.flip();
      remesher.smooth();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "# Iteration " << it_nb << " done : "
        << remesher.triangulation().number_of_vertices()
        << " vertices #" << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
      std::ostringstream ossi;
      ossi << "statistics_" << it_nb << ".txt";
      Tetrahedral_remeshing::internal::compute_statistics(
          remesher.triangulation(),
          remesher.imaginary_index(), cell_select, ossi.str().c_str());
#endif
    }

    std::size_t nb_extra_iterations = 3;
    while (it_nb++ < max_it + nb_extra_iterations)
    {
//      remesher.flip();
//      remesher.smooth();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "# Iteration " << it_nb << " (flip and smooth only) done : "
        << remesher.triangulation().number_of_vertices()
        << " vertices #" << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
      std::ostringstream ossi;
      ossi << "statistics_" << it_nb << ".txt";
      Tetrahedral_remeshing::internal::compute_statistics(
        remesher.triangulation(),
        remesher.imaginary_index(), cell_select, ossi.str().c_str());
#endif
    }

    remesher.postprocess();

    remesher.finalize();
    //remesher.triangulation() is now empty

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    const double angle_bound = 5.0;
    Tetrahedral_remeshing::debug::dump_cells_with_small_dihedral_angle(tr,
      angle_bound, remesher.imaginary_index(), cell_select, "bad_cells.mesh");
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    Tetrahedral_remeshing::internal::compute_statistics(tr,
      remesher.imaginary_index(), cell_select, "statistics_end.txt");
#endif
  }


  //template<typename Triangulation, typename NamedParameters>
  //void tetrahedral_adaptive_remeshing(Triangulation& tr,
  //                                    const double& target_edge_length,
  //                                    const NamedParameters& np)
  //{
  //  typedef typename Triangulation::Geom_traits K;
  //  CGAL::Uniform_sizing_field<K> sizing_field(target_edge_length);
  //  tetrahedral_adaptive_remeshing(tr, sizing_field, np);
  //}

  //template<typename Triangulation, typename SizingField>
  //void tetrahedral_adaptive_remeshing(Triangulation& tr,
  //                                    const SizingField& sizing_field)
  //{
  //  tetrahedral_adaptive_remeshing(tr, sizing_field,
  //    Polygon_mesh_processing::parameters::all_default());
  //}

  template<typename Triangulation>
  void tetrahedral_adaptive_remeshing(Triangulation& tr,
                                      const double& target_edge_length)
  {
    tetrahedral_adaptive_remeshing(tr, target_edge_length,
      Polygon_mesh_processing::parameters::all_default());
  }

}//end namespace CGAL

#endif //TETRAHEDRAL_REMESHING_H
