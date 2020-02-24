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

#include <CGAL/Mesh_3/config.h>
#include <boost/parameter/preprocessor.hpp>
#include <CGAL/boost/parameter.h>
#include <CGAL/Mesh_edge_criteria_3.h>
#include <CGAL/Mesh_facet_criteria_3.h>
#include <CGAL/Mesh_cell_criteria_3.h>
#include <cfloat> // for the macro DBL_MAX
#include <boost/type_traits/is_base_of.hpp>
namespace CGAL {
  
namespace parameters {

// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/boost/parameter.h>
CGAL_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

  BOOST_PARAMETER_NAME( (edge_size, tag) edge_size_ )
  BOOST_PARAMETER_NAME( (edge_sizing_field, tag) edge_sizing_field_ )
  BOOST_PARAMETER_NAME( (facet_angle, tag) facet_angle_ )
  BOOST_PARAMETER_NAME( (facet_size, tag) facet_size_ )
  BOOST_PARAMETER_NAME( (facet_sizing_field, tag) facet_sizing_field_ )
  BOOST_PARAMETER_NAME( (facet_distance, tag) facet_distance_ )
  BOOST_PARAMETER_NAME( (facet_topology, tag) facet_topology_ )
  BOOST_PARAMETER_NAME( (cell_radius_edge, tag) cell_radius_edge_ )
  BOOST_PARAMETER_NAME( (cell_radius_edge_ratio, tag) cell_radius_edge_ratio_ )
  BOOST_PARAMETER_NAME( (cell_size, tag) cell_size_ )
  BOOST_PARAMETER_NAME( (cell_sizing_field, tag) cell_sizing_field_ )
  BOOST_PARAMETER_NAME( (sizing_field, tag) sizing_field_ )

CGAL_PRAGMA_DIAG_POP

} // end namespace parameters
  
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
  template <class ArgumentPack>
  Mesh_criteria_3_impl(const ArgumentPack& args)
    : edge_criteria_(args[parameters::edge_size
                          | args[parameters::edge_sizing_field
                                 | args[parameters::sizing_field | FT(DBL_MAX)] ] ])
    , facet_criteria_(args[parameters::facet_angle | FT(0)],
                      args[parameters::facet_size
                           | args[parameters::facet_sizing_field
                                  | args[parameters::sizing_field | FT(0)] ] ],
                      args[parameters::facet_distance | FT(0)],
                      args[parameters::facet_topology | CGAL::FACET_VERTICES_ON_SURFACE])
    , cell_criteria_(args[parameters::cell_radius_edge_ratio
                          | args[parameters::cell_radius_edge | FT(0)] ],
                     args[parameters::cell_size
                          | args[parameters::cell_sizing_field
                                 | args[parameters::sizing_field | FT(0)] ] ])
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
  
  // For convenient constructor call (see examples)
  BOOST_PARAMETER_CONSTRUCTOR(Mesh_criteria_3, (Base), parameters::tag,
                              (optional (edge_size_,*)
                                        (edge_sizing_field_,*)
                                        (facet_angle_,*)
                                        (facet_size_,*)
                                        (facet_sizing_field_,*)
                                        (facet_distance_,*)
                                        (facet_topology_,*)
                                        (cell_radius_edge_,*)
                                        (cell_size_,*)
                                        (cell_sizing_field_,*)
                                        (sizing_field_,*)
                              ))
  
};  // end class Mesh_criteria_3

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_CRITERIA_3_H
