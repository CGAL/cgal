// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef TETRAHEDRAL_REMESHING_H
#define TETRAHEDRAL_REMESHING_H

#include <CGAL/license/Tetrahedral_remeshing.h>

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
* It is recommended to use `CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3`
* for the first parameter, as it inherits from `Triangulation_3` with a `TDS` suitable for this function.
*
* This function takes as input a 3-dimensional triangulation
* and performs a sequence of atomic operations
* in order to generate as output a high quality mesh with a prescribed
* uniform density.
* These atomic operations are performed as follows:
*   - edge splits, until all edges satisfy a prescribed length criterion,
*   - edge collapses, until all edges satisfy a prescribed length criterion,
*   - edge flips, to locally improve dihedral angles, until they can
*     no longer be improved by flipping,
*   - global smoothing by vertex relocations,
*   - re-projection of boundary vertices to the initial surface.
*
* This remeshing function can deal with multi-domains, boundaries, and features.
* It preserves the geometry of
* subdomains throughout the remeshing process.
*
* Subdomains are defined by indices that
* are stored in the cells of the input triangulation, following the `MeshCellBase_3`
* concept (refined by `RemeshingCellBase_3`).
* The surfacic interfaces between subdomains are formed by facets whose two incident cells
* have different subdomain indices.
* The edges where three or more subdomains meet form feature polylines,
* and are considered as constrained edges.
*
*
* @tparam Traits is the geometric traits, model of `RemeshingTriangulationTraits_3`
* @tparam TDS is the triangulation data structure for `Triangulation_3`,
*             model of ` TriangulationDataStructure_3`,
*             with cell base model of `RemeshingCellBase_3`
*             and vertex base model of `RemeshingVertexBase_3`.
* @tparam SLDS is an optional parameter for `Triangulation_3`, that
*             specifies the type of the spatial lock data structure.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param tr the triangulation to be remeshed, of type `Triangulation_3<Traits, TDS, SLDS>`.
*           `Remeshing_triangulation` is a helper class that satisfies all the requirements
*           of its template parameters.
* @param target_edge_length the uniform target edge length. This parameter provides a
*          mesh density target for the remeshing algorithm.
* @param np optional sequence of \ref bgl_namedparameters "Named Parameters"
*          among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{the number of iterations for the full sequence of atomic operations
*                           performed (listed in the above description)}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`1`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{remesh_boundaries}
*     \cgalParamDescription{If `false`, none of the input volume boundaries can be modified.
*                           Otherwise, the topology is preserved, but atomic operations
*                           can be performed on the surfaces, and along feature polylines,
*                           such that boundaries are remeshed.}
*     \cgalParamType{`bool`}
*     \cgalParamDefault{`true`}
*     \cgalParamExtra{Boundaries are between the exterior and the interior,
*                     between two subdomains, between the areas selected or not for remeshing
*                     (cf `Remeshing_cell_is_selected_map`), or defined
*                     by `Remeshing_edge_is_constrained_map` and `Remeshing_facet_is_constrained_map`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{edge_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tr`.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `std::pair<Triangulation_3::Vertex_handle, Triangulation_3::Vertex_handle>`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no edge is constrained}
*     \cgalParamExtra{A constrained edge can be split or collapsed, but not flipped.}
*     \cgalParamExtra{The pairs must be ordered to ensure consistency.}
*     \cgalParamExtra{During the meshing process, the set of constrained edges evolves consistently
*                     with edge splits and collapses, so the property map must be writable.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{facet_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each facet of `tr`.}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `Triangulation_3::Facet`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no facet is constrained}
*     \cgalParamExtra{A constrained facet can be split or collapsed, but not flipped.}
*     \cgalParamExtra{This map, contrary to the others, is not updated throughout the remeshing process.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{cell_is_selected_map}
*     \cgalParamDescription{a property map containing the selected - or - not status for each cell of `tr` for remeshing.
*                           Only selected cells are modified (and possibly their neighbors if surfaces are modified) by remeshing.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `Triangulation_3::Cell_handle`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where all cells of the domain
*                       (i.e. with a non-zero `Subdomain_index`) are selected.}
*     \cgalParamExtra{During the meshing process, the set of selected cells evolves consistently with
*                     the atomic operations that are performed, so the property map must be writable.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{smooth_constrained_edges}
*     \cgalParamDescription{If `true`, the end vertices of the edges set as
*                           constrained in `edge_is_constrained_map` move along the
*                           constrained polylines they belong to.}
*     \cgalParamType{`bool`}
*     \cgalParamDefault{`false`}
*     \cgalParamExtra{The endvertices of constraints listed
*                     by `edge_is_constrained_map`, and edges incident to at least three subdomains
*                     are made eligible to one dimensional smoothing, along the constrained polylines they belong to.
*                     Corners (i.e. vertices incident to more than 2 constrained edges) are not allowed
*                     to move at all.\n
*                     Note that activating the smoothing step on polyline constraints tends to reduce
*                     the quality of the minimal dihedral angle in the mesh.\n
*                     If `remesh_boundaries` is set to `false`, this parameter is ignored.}
*   \cgalParamNEnd
*
*  \cgalNamedParamsEnd
*
* \sa `CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3`
*
* @todo implement non-uniform sizing field instead of uniform target edge length
*/
template<typename Traits, typename TDS, typename SLDS,
         typename NamedParameters>
void tetrahedral_isotropic_remeshing(
  CGAL::Triangulation_3<Traits, TDS, SLDS>& tr,
  const double& target_edge_length,
  const NamedParameters& np)
{
  typedef CGAL::Triangulation_3<Traits, TDS, SLDS> Triangulation;
  tetrahedral_isotropic_remeshing(
    tr,
    [target_edge_length](const typename Triangulation::Point& /* p */)
                        {return target_edge_length;},
    np);
}

template<typename Traits, typename TDS, typename SLDS,
         typename NamedParameters>
void tetrahedral_isotropic_remeshing(
  CGAL::Triangulation_3<Traits, TDS, SLDS>& tr,
  const float& target_edge_length,
  const NamedParameters& np)
{
  typedef CGAL::Triangulation_3<Traits, TDS, SLDS> Triangulation;
  tetrahedral_isotropic_remeshing(
    tr,
    [target_edge_length](const typename Triangulation::Point& /* p */)
                        {return target_edge_length; },
    np);
}

template<typename Traits, typename TDS, typename SLDS,
         typename SizingFunction,
         typename NamedParameters>
void tetrahedral_isotropic_remeshing(
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
  bool smooth_constrained_edges
    = choose_parameter(get_parameter(np, internal_np::smooth_constrained_edges),
                       false);

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
    internal_np::visitor_t,
    NamedParameters,
    Tetrahedral_remeshing::internal::Default_remeshing_visitor
  > ::type Visitor;
  Visitor visitor
    = choose_parameter(get_parameter(np, internal_np::visitor),
                       Tetrahedral_remeshing::internal::Default_remeshing_visitor());

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  std::cout << "Tetrahedral remeshing ("
            << "nb_iter = " << max_it << ", "
            << "protect = " << std::boolalpha << protect
            << ")" << std::endl;

  std::cout << "Init tetrahedral remeshing...";
  std::cout.flush();
#endif

  typedef Tetrahedral_remeshing::internal::Adaptive_remesher<
    Tr, SizingFunction, ECMap, FCMap, SelectionFunctor, Visitor> Remesher;
  Remesher remesher(tr, sizing, protect
                  , ecmap, fcmap
                  , smooth_constrained_edges
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

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  const double angle_bound = 5.0;
  Tetrahedral_remeshing::debug::dump_cells_with_small_dihedral_angle(tr,
    angle_bound, cell_select, "bad_cells.mesh");
  Tetrahedral_remeshing::internal::compute_statistics(tr,
    cell_select, "statistics_end.txt");
#endif
}

template<typename Traits, typename TDS, typename SLDS>
void tetrahedral_isotropic_remeshing(
  CGAL::Triangulation_3<Traits, TDS, SLDS>& tr,
  const double& target_edge_length)
{
  tetrahedral_isotropic_remeshing(tr, target_edge_length,
                                 CGAL::parameters::all_default());
}

/*!
* \ingroup PkgTetrahedralRemeshingRef
* converts the triangulation contained in the input to a `Triangulation_3`.
*
* This function should be used to generate a valid triangulation
* for tetrahedral remeshing, when the input triangulation is generated with the
* tetrahedral mesh generation package.
*
* @tparam Tr is the underlying triangulation for `Mesh_complex_3_in_triangulation_3`.
*            It can be instantiated with any 3D regular triangulation of CGAL provided
*            that its vertex and cell base classes are models of the concepts
*            `MeshVertexBase_3` (refined by `RemeshingCellBase_3`)
*            and `MeshCellBase_3` (refined by `RemeshingVertexBase_3`), respectively.
* @tparam CornerIndex is the type of the indices for feature corners.
*            If `c3t3` has been generated using `CGAL::make_mesh_3()`, it must match
*            the `Corner_index` type of the model of the `MeshDomainWithFeatures_3` concept used for mesh generation.
* @tparam CurveIndex is the type of the indices for feature curves.
*            If `c3t3` has been generated using `CGAL::make_mesh_3()`, it must match
*            the `Curve_index` type of the model of the `MeshDomainWithFeatures_3` concept used for mesh generation.
*
* @param c3t3 the complex containing the triangulation to be remeshed.
*/

template<typename Tr,
         typename CornerIndex,
         typename CurveIndex>
CGAL::Triangulation_3<typename Tr::Geom_traits,
                      typename Tr::Triangulation_data_structure>
convert_to_triangulation_3(
  CGAL::Mesh_complex_3_in_triangulation_3<Tr, CornerIndex, CurveIndex> c3t3)
{
  using GT   = typename Tr::Geom_traits;
  using TDS  = typename Tr::Triangulation_data_structure;

  CGAL::Triangulation_3<GT, TDS> tr;
  tr.swap(c3t3.triangulation());
  return tr;
}

///////////////////////////////////////////////////
/////// MESH_COMPLEX_3_IN_TRIANGULATION_3 /////////
///////////////////////////////////////////////////

///////
////// Warning with using this version :
////// the triangulation after remeshing is not regular anymore
////// the empty-sphere property is not respected
///////
template<typename Tr,
         typename CornerIndex, typename CurveIndex,
         typename NamedParameters>
void tetrahedral_isotropic_remeshing(
  CGAL::Mesh_complex_3_in_triangulation_3<Tr, CornerIndex, CurveIndex>& c3t3,
  const double& target_edge_length,
  const NamedParameters& np)
{
  tetrahedral_isotropic_remeshing(
    c3t3,
    [target_edge_length](const typename Tr::Point& /* p */)
                        {return target_edge_length; },
    np);
}

template<typename Tr,
         typename CornerIndex, typename CurveIndex,
         typename NamedParameters>
void tetrahedral_isotropic_remeshing(
  CGAL::Mesh_complex_3_in_triangulation_3<Tr, CornerIndex, CurveIndex>& c3t3,
  const float& target_edge_length,
  const NamedParameters& np)
{
  tetrahedral_isotropic_remeshing(
    c3t3,
    [target_edge_length](const typename Tr::Point& /*p*/)
                        {return target_edge_length; },
    np);
}

template<typename Tr,
         typename CornerIndex, typename CurveIndex>
void tetrahedral_isotropic_remeshing(
  CGAL::Mesh_complex_3_in_triangulation_3<Tr, CornerIndex, CurveIndex>& c3t3,
  const double& target_edge_length)
{
  return tetrahedral_isotropic_remeshing(c3t3, target_edge_length,
                                        CGAL::parameters::all_default());
}

template<typename Tr,
         typename CornerIndex, typename CurveIndex,
         typename SizingFunction,
         typename NamedParameters>
void tetrahedral_isotropic_remeshing(
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

  bool smooth_constrained_edges
    = choose_parameter(get_parameter(np, internal_np::smooth_constrained_edges),
      false);

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
              internal_np::visitor_t,
              NamedParameters,
              Tetrahedral_remeshing::internal::Default_remeshing_visitor
              > ::type Visitor;
  Visitor visitor
    = choose_parameter(get_parameter(np, internal_np::visitor),
                       Tetrahedral_remeshing::internal::Default_remeshing_visitor());

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  std::cout << "Tetrahedral remeshing ("
            << "nb_iter = " << max_it << ", "
            << "protect = " << std::boolalpha << protect
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
                    , smooth_constrained_edges
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

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
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
