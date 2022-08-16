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
#include <boost/type_traits/is_base_of.hpp>
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
    :edge_criteria_(parameters::choose_parameter(parameters::get_parameter(np, internal_np::edge_size_param),
                                                 parameters::choose_parameter(parameters::get_parameter(np, internal_np::edge_sizing_field_param),
                                                                              parameters::choose_parameter(parameters::get_parameter(np, internal_np::sizing_field_param), FT(DBL_MAX))))),
    facet_criteria_(parameters::choose_parameter(parameters::get_parameter(np, internal_np::facet_angle_param), FT(0)),
                    parameters::choose_parameter(parameters::get_parameter(np, internal_np::facet_size_param),
                                                 parameters::choose_parameter(parameters::get_parameter(np, internal_np::facet_sizing_field_param),
                                                                              parameters::choose_parameter(parameters::get_parameter(np, internal_np::sizing_field_param), FT(0)))),
                    parameters::choose_parameter(parameters::get_parameter(np, internal_np::facet_distance_param), FT(0)),
                    parameters::choose_parameter(parameters::get_parameter(np, internal_np::facet_topology_param), CGAL::FACET_VERTICES_ON_SURFACE)),
    cell_criteria_(parameters::choose_parameter(parameters::get_parameter(np, internal_np::cell_radius_edge_ratio_param),
                                                parameters::choose_parameter(parameters::get_parameter(np, internal_np::cell_radius_edge_param), FT(0))),
                   parameters::choose_parameter(parameters::get_parameter(np, internal_np::cell_size_param),
                                                parameters::choose_parameter(parameters::get_parameter(np, internal_np::cell_sizing_field_param),
                                                                             parameters::choose_parameter(parameters::get_parameter(np, internal_np::sizing_field_param), FT(0)))))
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
    CGAL_static_assertion((boost::is_base_of<
                           typename Facet_criteria::Abstract_criterion,
                           Facet_criterion
                           >::value));
    facet_criteria_.add(criterion);
  }

  template <typename Cell_criterion>
  void add_cell_criterion(Cell_criterion* criterion) {
    CGAL_static_assertion((boost::is_base_of<
                           typename Cell_criteria::Abstract_criterion,
                           Cell_criterion
                           >::value));
    cell_criteria_.add(criterion);
  }

private:
  Edge_criteria edge_criteria_;
  Facet_criteria facet_criteria_;
  Cell_criteria cell_criteria_;

};  // end class Mesh_criteria_3_impl

} // end namespace internal



// Class Mesh_criteria_3
// Provides default mesh criteria to drive Mesh_3 process
template <typename Tr,
          typename EdgeCriteria = Mesh_edge_criteria_3<Tr>,
          typename FacetCriteria = Mesh_facet_criteria_3<Tr>,
          typename CellCriteria = Mesh_cell_criteria_3<Tr> >
class Mesh_criteria_3
  : public internal::Mesh_criteria_3_impl< Tr,
                                           EdgeCriteria,
                                           FacetCriteria,
                                           CellCriteria >
{
  typedef internal::Mesh_criteria_3_impl< Tr,
                                          EdgeCriteria,
                                          FacetCriteria,
                                          CellCriteria>   Base;

public:
  typedef typename Base::Edge_criteria    Edge_criteria;
  typedef typename Base::Facet_criteria   Facet_criteria;
  typedef typename Base::Cell_criteria    Cell_criteria;

  // Constructor
  Mesh_criteria_3(Facet_criteria facet_criteria,
                  Cell_criteria cell_criteria)
    : Base(facet_criteria,
           cell_criteria) {}

  // Constructor
  Mesh_criteria_3(Edge_criteria edge_criteria,
                  Facet_criteria facet_criteria,
                  Cell_criteria cell_criteria)
    : Base(edge_criteria,
           facet_criteria,
           cell_criteria) {}
/*!
      \brief Construction from criteria parameters.
      This constructor uses named
      parameters (from <I>Boost.Parameter</I>) for convenient criteria
      construction.

         \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

         \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below:

         The following are optional named parameters.

         \cgalNamedParamsBegin
           \cgalParamNBegin{edge_size}
             \cgalParamDescription{a scalar field (resp. a constant) providing a space varying
                                  (resp. a uniform)
                                  upper bound for the lengths of curve edges. This parameter has to be set to a positive
                                  value when 1-dimensional features protection is used.}

           \cgalParamNBegin{facet_angle}
             \cgalParamDescription{a lower bound for the angles (in degrees) of the
                                   surface mesh facets.}

           \cgalParamNBegin{facet_size}
             \cgalParamDescription{ a scalar field (resp. a constant) describing
                                    a space varying (resp. a uniform) upper-bound or for the radii of the surface Delaunay balls.}

           \cgalParamNBegin{facet_distance}
             \cgalParamDescription{ a scalar field (resp. a constant) describing a space varying (resp. a uniform)
                                    upper bound for the distance between the facet circumcenter and the center of its surface
                                    Delaunay ball.}
           \cgalParamNBegin{facet_topology}
             \cgalParamDescription{ the set of topological constraints
                                    which have to be verified by each surface facet. See `Mesh_facet_topology` manual page to
                                    get all possible values.}
              \cgalParamDefault{CGAL::FACET_VERTICES_ON_SURFACE}

            \cgalParamNBegin{cell_radius_edge_ratio}
             \cgalParamDescription{ an upper bound for the radius-edge ratio of the mesh tetrahedra.}

             \cgalParamNBegin{cell_size}
             \cgalParamDescription{ a scalar field (resp. a constant) describing
                                    a space varying (resp. a uniform) upper-bound for the circumradii of the mesh tetrahedra.}


         \cgalNamedParamsEnd
       */

template<typename CGAL_NP_TEMPLATE_PARAMETERS>
Mesh_criteria_3(const CGAL_NP_CLASS& np = parameters::default_values()): Base(np)
{
}


template<typename ... CGAL_NP_TEMPLATE_PARAMETERS_VARIADIC>
Mesh_criteria_3(const CGAL_NP_CLASS& ... nps):Mesh_criteria_3(internal_np::combine_named_parameters(nps...))
{
}

};  // end class Mesh_criteria_3

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_CRITERIA_3_H
