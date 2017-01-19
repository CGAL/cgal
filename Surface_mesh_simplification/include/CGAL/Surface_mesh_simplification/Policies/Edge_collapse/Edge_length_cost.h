// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
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
  template<class ECM>
class Edge_length_cost
{
public:
  /*  
  typedef ECM_ ECM ;
  
  typedef Edge_profile<ECM> Profile ;
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
 
