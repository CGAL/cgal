// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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

#include <CGAL/Mesh_3/global_parameters.h>
#include <CGAL/Mesh_facet_criteria_3.h>
#include <CGAL/Mesh_cell_criteria_3.h>

namespace CGAL {
  
namespace parameters {
  BOOST_PARAMETER_NAME( (facet_angle, tag) facet_angle_ )
  BOOST_PARAMETER_NAME( (facet_size, tag) facet_size_ )
  BOOST_PARAMETER_NAME( (facet_distance, tag) facet_distance_ )
  BOOST_PARAMETER_NAME( (cell_radius_edge, tag) cell_radius_edge_ )
  BOOST_PARAMETER_NAME( (cell_size, tag) cell_size_ )
} // end namespace parameters
  
namespace internal {

// Class Mesh_criteria_3_impl
template <class Tr, typename FacetCriteria, typename CellCriteria>
class Mesh_criteria_3_impl
{
public:
  typedef FacetCriteria     Facet_criteria;
  typedef CellCriteria      Cell_criteria;
  
  // Constructor
  Mesh_criteria_3_impl(const Facet_criteria& facet_criteria,
                       const Cell_criteria& cell_criteria)
    : facet_criteria_(facet_criteria)
    , cell_criteria_(cell_criteria) {}
  
  // This template constructor is not instantiated when named parameters
  // are not used, so Facet_criteria and Cell_criteria construction from FT
  // is not a problem
  template <class ArgumentPack>
  Mesh_criteria_3_impl(const ArgumentPack& args)
    : facet_criteria_(args[parameters::facet_angle | 0],
                      args[parameters::facet_size | 0],
                      args[parameters::facet_distance | 0])
    , cell_criteria_(args[parameters::cell_radius_edge | 0],
                     args[parameters::cell_size | 0] )              { }
  
  const Facet_criteria& facet_criteria() const { return facet_criteria_; }
  const Cell_criteria& cell_criteria() const { return cell_criteria_; }
  
private:
  Facet_criteria facet_criteria_;
  Cell_criteria cell_criteria_;
  
};  // end class Mesh_criteria_3_impl  

} // end namespace internal
  
  
  
// Class Mesh_criteria_3
// Provides default mesh criteria to drive Mesh_3 process
template <class Tr,
          typename FacetCriteria = Mesh_facet_criteria_3<Tr>,
          typename CellCriteria = Mesh_cell_criteria_3<Tr> >
class Mesh_criteria_3
  : public internal::Mesh_criteria_3_impl<Tr, FacetCriteria, CellCriteria>
{
  typedef internal::Mesh_criteria_3_impl<Tr, FacetCriteria, CellCriteria> Base;
  
public:
  typedef typename Base::Facet_criteria   Facet_criteria;
  typedef typename Base::Cell_criteria    Cell_criteria;
  
  // Constructor
  Mesh_criteria_3(const Facet_criteria& facet_criteria,
                  const Cell_criteria& cell_criteria)
    : Base(facet_criteria, cell_criteria) {}
  
  // For convenient constructor call (see examples)
  BOOST_PARAMETER_CONSTRUCTOR(Mesh_criteria_3, (Base), parameters::tag,
                              (optional (facet_angle_,*)
                                        (facet_size_,*)
                                        (facet_distance_,*)
                                        (cell_radius_edge_,*)
                                        (cell_size_,*)           ))
  
};  // end class Mesh_criteria_3

}  // end namespace CGAL


#endif // CGAL_MESH_CRITERIA_3_H
