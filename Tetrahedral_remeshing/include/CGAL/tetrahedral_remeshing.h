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
  * remeshes a tetrahedral mesh.
  *
  * This operation sequentially performs edge splits, edge collapses,
  * edge flips, smoothing and projection to the initial surface to generate
  * a quality mesh with a prescribed edge length.
  *
  * @tparam Triangulation model of `Triangulation_3`,
  * with cell base model of `TriangulationCellBaseWithInfo_3`
  * 
  * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"

  * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
  * \cgalNamedParamsBegin
  *  \cgalParamBegin{protect_boundaries} If `true`, the
  *     volume boundaries cannot be modified (no modification of boundaries in this version)
  *  \cgalParamEnd
  *  \cgalParamBegin{number_of_iterations} the number of iterations for the sequence of atomic operations
  *     performed (listed in the above description)
  *  \cgalParam
  *  \cgalParamBegin{cell_selector} a functor that returns a boolean setting whether the given
  *    `Triangulation::Cell_handle` should be part of the remeshing (by default, cells are all part
  *    of the remeshing)
  *  \cgalParamEnd
  * \cgalNamedParamsEnd
  */

  //  * @tparam SizingField model of `CGAL::Sizing_field`

  //*  \cgalParamBegin{ adaptive } If `true`, size of elements adapts
  //*     ....
  //*  \cgalParamEnd
    //*  \cgalParamBegin{ edge_is_constrained_map } a property map containing the
    //*    constrained - or - not status of each edge of `tr`. A constrained edge can be split
    //*    or collapsed, but not flipped, nor its endpoints moved by smoothing.
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
    CGAL_assertion(tr.is_valid(true));

    typedef Triangulation Tr;
    typedef typename Tr::Edge Edge;

    using boost::choose_param;
    using boost::get_param;

    bool protect = choose_param(get_param(np, internal_np::protect_boundaries),
                                true);
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

    typedef Tetrahedral_remeshing::internal::No_constraint_pmap<Edge> No_constraint;

    typedef typename boost::lookup_named_param_def <
      internal_np::edge_is_constrained_t,
      NamedParameters,
      No_constraint//default
    > ::type ECMap;
    ECMap ecmap = choose_param(get_param(np, internal_np::edge_is_constrained)
                             , No_constraint());

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Init tetrahedral remeshing...";
    std::cout.flush();
#endif

    typedef Tetrahedral_remeshing::internal::Adaptive_remesher<
      Tr, ECMap, SelectionFunctor> Remesher;
    Remesher remesher(tr, target_edge_length, protect, ecmap
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
    while (it_nb++ < max_it && !remesher.resolution_reached())
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "# Iteration " << it_nb << " #" << std::endl;
#endif
      remesher.split();
      remesher.collapse();
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
      remesher.flip();
      remesher.smooth();

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
