// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H

#include <CGAL/license/Surface_mesh_simplification.h>


#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

namespace CGAL {

namespace Surface_mesh_simplification
{

//
// Edge-length cost: the squared length of the collapsing edge
//
  template<class TM>
class Edge_length_cost
{
public:
  /*  
  typedef TM_ TM ;
  
  typedef Edge_profile<TM> Profile ;
  typedef typename Profile::Point Point;  
  typedef typename Kernel_traits<Point>::Kernel Kernel ;
  typedef typename Kernel::FT FT ;
  typedef optional<FT> result_type ;
  */
public:

  Edge_length_cost()
  {}

  template <typename Profile, typename T> 
  optional<typename Profile::FT> operator()( Profile const& aProfile, T const& /*aPlacement*/ ) const
  {
    typedef optional<typename Profile::FT> result_type;
    return result_type(squared_distance(aProfile.p0(),aProfile.p1()));
  }
  
};


} // namespace Surface_mesh_simplification


} //namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H
// EOF //
 
