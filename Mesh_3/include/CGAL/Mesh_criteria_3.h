// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Implements default meshing criteria to drive Mesh_3 process
//******************************************************************************


#ifndef CGAL_MESH_CRITERIA_3_H
#define CGAL_MESH_CRITERIA_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Mesh_3/config.h>
#include <CGAL/Mesh_edge_criteria_3.h>
#include <CGAL/Mesh_facet_criteria_3.h>
#include <CGAL/Mesh_cell_criteria_3.h>
#include <cfloat> // for the macro DBL_MAX

namespace CGAL {


namespace internal {

// Class Mesh_criteria_3_impl
template < typename Tr,
           typename EdgeCriteria,
           typename FacetCriteria,
           typename CellCriteria >
class Mesh_criteria_3_impl
{
  typedef typename Tr::Geom_traits::FT FT;

public:
  typedef EdgeCriteria      Edge_criteria;
  typedef FacetCriteria     Facet_criteria;
  typedef CellCriteria      Cell_criteria;

  // Constructor
  Mesh_criteria_3_impl(Facet_criteria facet_criteria,
                       Cell_criteria cell_criteria)
    : edge_criteria_(0)
    , facet_criteria_(facet_criteria)
    , cell_criteria_(cell_criteria)
  { }

  // Constructor
  Mesh_criteria_3_impl(Edge_criteria edge_criteria,
                       Facet_criteria facet_criteria,
                       Cell_criteria cell_criteria)
    : edge_criteria_(edge_criteria)
    , facet_criteria_(facet_criteria)
    , cell_criteria_(cell_criteria)
  { }

  // This template constructor is not instantiated when named parameters
  // are not used, so Facet_criteria and Cell_criteria construction from FT
  // is not a problem
  template <typename CGAL_NP_TEMPLATE_PARAMETERS>
  Mesh_criteria_3_impl(const CGAL_NP_CLASS& np)
    :edge_criteria_(parameters::choose_parameter(parameters::get_parameter_reference(np, internal_np::edge_size_param), FT(DBL_MAX)),
                    parameters::choose_parameter(parameters::get_parameter(np, internal_np::edge_min_size_param), FT(0)),
                    parameters::choose_parameter(parameters::get_parameter_reference(np, internal_np::edge_distance_param), FT(0))),
    facet_criteria_(parameters::choose_parameter(parameters::get_parameter(np, internal_np::facet_angle_param), FT(0)),
                    parameters::choose_parameter(parameters::get_parameter_reference(np, internal_np::facet_size_param), FT(0)),
                    parameters::choose_parameter(parameters::get_parameter_reference(np, internal_np::facet_distance_param), FT(0)),
                    parameters::choose_parameter(parameters::get_parameter(np, internal_np::facet_topology_param), CGAL::FACET_VERTICES_ON_SURFACE),
                    parameters::choose_parameter(parameters::get_parameter(np, internal_np::facet_min_size_param), FT(0))),
    cell_criteria_(parameters::choose_parameter(parameters::get_parameter(np, internal_np::cell_radius_edge_ratio_param), FT(0)),
                   parameters::choose_parameter(parameters::get_parameter_reference(np, internal_np::cell_size_param), FT(0)),
                   parameters::choose_parameter(parameters::get_parameter(np, internal_np::cell_min_size_param), FT(0)))
  { }

#ifndef CGAL_NO_DEPRECATED_CODE
  const Edge_criteria& edge_criteria() const { return edge_criteria_; }
  const Facet_criteria& facet_criteria() const { return facet_criteria_; }
  const Cell_criteria& cell_criteria() const { return cell_criteria_; }
#endif

  const Edge_criteria& edge_criteria_object() const { return edge_criteria_; }
  const Facet_criteria& facet_criteria_object() const { return facet_criteria_; }
  const Cell_criteria& cell_criteria_object() const { return cell_criteria_; }

  template <typename Facet_criterion>
  void add_facet_criterion(Facet_criterion* criterion) {
    static_assert(std::is_base_of<
                           typename Facet_criteria::Abstract_criterion,
                           Facet_criterion
                           >::value);
    facet_criteria_.add(criterion);
  }

  template <typename Cell_criterion>
  void add_cell_criterion(Cell_criterion* criterion) {
    static_assert(std::is_base_of<
                           typename Cell_criteria::Abstract_criterion,
                           Cell_criterion
                           >::value);
    cell_criteria_.add(criterion);
  }

private:
  Edge_criteria edge_criteria_;
  Facet_criteria facet_criteria_;
  Cell_criteria cell_criteria_;

};  // end class Mesh_criteria_3_impl

} // end namespace internal




/*!
\ingroup PkgMesh3MeshClasses

The class gathers the refinement criteria for mesh tetrahedra and
surface facets where
surface facets are facets in the mesh approximating the domain surface patches.
In addition, for domains with exposed 1-dimensional features,
the class `Mesh_criteria_3`
handles the definition of a sizing field to guide the discretization of
1-dimensional features.

\tparam Tr has to be instantiated with the type used for
`C3T3::Triangulation`,
where `C3T3` is the model of `MeshComplex_3InTriangulation_3`
used in the mesh generation process,
and `C3T3::Triangulation` its nested triangulation type.

\cgalModels{MeshCriteria_3,MeshCriteriaWithFeatures_3}

\cgalHeading{Example}

\code{.cpp}

// Create a Mesh_criteria_3<Tr> object with all cell and facet parameters set
Mesh_criteria_3<Tr> criteria (parameters::facet_angle(30).
                              parameters::facet_size(1).
                              parameters::facet_distance(0.1).
                              parameters::cell_radius_edge_ratio(2).
                              parameters::cell_size(1.5));

// Create a Mesh_criteria_3<Tr> object with size ignored (note that the order changed)
Mesh_criteria_3<Tr> criteria (parameters::cell_radius_edge_ratio(2).
                              parameters::facet_angle(30).
                              parameters::facet_distance(0.1));

\endcode

\sa `MeshCriteria_3`
\sa `MeshCriteriaWithFeatures_3`
\sa `MeshCellCriteria_3`
\sa `MeshEdgeCriteria_3`
\sa `MeshFacetCriteria_3`
\sa `MeshDomainField_3`
\sa `CGAL::Mesh_cell_criteria_3<Tr>`
\sa `CGAL::Mesh_edge_criteria_3<Tr>`
\sa `CGAL::Mesh_facet_criteria_3<Tr>`
\sa `CGAL::Mesh_facet_topology`

*/
template <typename Tr,
          typename EdgeCriteria = Mesh_edge_criteria_3<Tr>,
          typename FacetCriteria = Mesh_facet_criteria_3<Tr>,
          typename CellCriteria = Mesh_cell_criteria_3<Tr> >
class Mesh_criteria_3
#ifndef DOXYGEN_RUNNING
  : public internal::Mesh_criteria_3_impl< Tr,
                                           EdgeCriteria,
                                           FacetCriteria,
                                           CellCriteria >
#endif
{
  typedef internal::Mesh_criteria_3_impl< Tr,
                                          EdgeCriteria,
                                          FacetCriteria,
                                          CellCriteria>   Base;

public:
#ifdef DOXYGEN_RUNNING
/// \name Types
/// @{

/*!
The criteria for edges.
*/
typedef Mesh_edge_criteria_3<Tr> Edge_criteria;

/*!
The criteria for facets.
*/
typedef Mesh_facet_criteria_3<Tr> Facet_criteria;

/*!
The
criteria for cells.
*/
typedef Mesh_cell_criteria_3<Tr> Cell_criteria;

/// @}
#else
  typedef typename Base::Edge_criteria    Edge_criteria;
  typedef typename Base::Facet_criteria   Facet_criteria;
  typedef typename Base::Cell_criteria    Cell_criteria;
#endif

  /// Construction from facet and cell criteria, the edge criteria are ignored in this case.
  Mesh_criteria_3(Facet_criteria facet_criteria,
                  Cell_criteria cell_criteria)
    : Base(facet_criteria,
           cell_criteria) {}

  /// Constructor from edge, face, and cell criteria.
  Mesh_criteria_3(Edge_criteria edge_criteria,
                  Facet_criteria facet_criteria,
                  Cell_criteria cell_criteria)
    : Base(edge_criteria,
           facet_criteria,
           cell_criteria) {}
/*!
 * \brief Construction from criteria parameters.
 *
 * Note that each size or distance parameter can be specified using two ways: either as a scalar field or
 * as a numerical value when the field is uniform.
 *
 * If not specified, each parameter has a default value such that the corresponding criterion will be ignored.
 * Numerical sizing or distance values, as well as scalar fields should be given in the unit used for coordinates
 * of points in the mesh domain class of the mesh generation process.
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{edge_size}
 *     \cgalParamDescription{a scalar field (resp. a constant) providing a space varying
 *                           (resp. a uniform)
 *                           upper bound for the lengths of curve edges. This parameter has to be set to a positive
 *                           value when 1-dimensional features protection is used.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{edge_min_size}
 *     \cgalParamDescription{a desired uniform lower-bound for the lengths of curve edges.
 *                           Only feature edges with a length larger than this bound will be refined.
 *                           If a feature edge is too small with respect to this criterion,
 *                           it will not be refined whether the other criteria are met or not.}
 *     \cgalParamExtra{If this criterion is applied during the meshing process,
 *                     the feature protection algorithm correctness is not guaranteed anymore,
 *                     and the output mesh may contain incorrect polyline features,
 *                     or have missing polyline features.}
 *     \cgalParamExtra{Note this lower-bound may not be respected everywhere in the output mesh,
 *                     to keep the feature graph valid.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{edge_distance}
 *     \cgalParamDescription{a scalar field (resp. a constant) describing a space varying
 *                           (resp. a uniform) upper bound for the distance between the
 *                           edge and its corresponding 1D feature.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{facet_angle}
 *     \cgalParamDescription{a lower bound for the angles (in degrees) of the
 *                            surface mesh facets.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{facet_size}
 *     \cgalParamDescription{a scalar field (resp. a constant) describing
 *                           a space varying (resp. a uniform) upper-bound or for the radii of the surface Delaunay balls.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{facet_min_size}
 *     \cgalParamDescription{a constant describing a uniform lower-bound for the radii of the surface Delaunay balls.
 *                           Only facets with a radius larger than this bound will be refined.
 *                           If a facet is too small with respect to this criterion,
 *                           it will not be refined whether the other criteria are met or not.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{facet_distance}
 *     \cgalParamDescription{a scalar field (resp. a constant) describing a space varying (resp. a uniform)
 *                           upper bound for the distance between the facet circumcenter and the center of its surface
 *                           Delaunay ball.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{facet_topology}
 *     \cgalParamDescription{the set of topological constraints
 *                           which have to be verified by each surface facet. See `Mesh_facet_topology` manual page to
 *                           get all possible values.}
 *     \cgalParamDefault{CGAL::FACET_VERTICES_ON_SURFACE}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{cell_radius_edge_ratio}
 *     \cgalParamDescription{ an upper bound for the radius-edge ratio of the mesh tetrahedra.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{cell_size}
 *     \cgalParamDescription{ a scalar field (resp. a constant) describing
 *                            a space varying (resp. a uniform) upper-bound for the circumradii of the mesh tetrahedra.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{cell_min_size}
 *     \cgalParamDescription{ a constant describing a uniform lower-bound for the radii of the circumradii
 *                            of the mesh tetrahedra.
 *                            Only tetrahedra with a circumradius larger than this bound will be refined.
 *                            If a cell is too small with respect to this criterion,
 *                            it will not be refined whether the other criteria are met or not.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 */
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_criteria_3(const CGAL_NP_CLASS& np = parameters::default_values()): Base(np)
{
}

// Overload handling parameters passed with operator=
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_1,
         typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT_2,
         typename ... NP>
Mesh_criteria_3(const CGAL_NP_CLASS_1&  np1,
                const CGAL_NP_CLASS_2&  np2,
                const NP& ... nps)
  : Mesh_criteria_3(internal_np::combine_named_parameters(np1, np2, nps...))
{
}

};  // end class Mesh_criteria_3

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_CRITERIA_3_H
