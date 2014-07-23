// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_MINIMUM_LENGTH_PREDICATE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_MINIMUM_LENGTH_PREDICATE_H 1

/// @cond CGAL_DOCUMENT_INTERNAL

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

namespace SMS = CGAL::Surface_mesh_simplification;

namespace CGAL {

namespace internal {

//*******************************************************************************************************************
//                                -= stopping condition predicate =-
//
// Determines whether the simplification has finished.
// The arguments are (current_cost,vertex,vertex,is_edge,initial_pair_count,current_pair_count,surface) and the result is bool
//
//*******************************************************************************************************************

// 
// Stops when all the edges have length greater than the minimum required edge length 
//
template<class ECM_>    
class Minimum_length_predicate
{
public:

  typedef ECM_ ECM ;

  typedef SMS::Edge_profile<ECM> Profile ;
  
private :

  typedef typename halfedge_graph_traits<ECM>::Point Point ;
  typedef typename Kernel_traits<Point>::Kernel      Kernel ;

public :
  
  typedef typename boost::graph_traits<ECM>::halfedge_descriptor halfedge_descriptor ;
  typedef typename boost::graph_traits<ECM>::edges_size_type size_type ;
  
  typedef typename Kernel::FT FT ;

public :
  
  Minimum_length_predicate( FT aThres ) : mThres(aThres * aThres) {}
  
  bool operator()( FT const&      aCurrentCost
                 , Profile const& //aEdgeProfile
                 , size_type      //aInitialCount
                 , size_type      //aCurrentCount
                 ) const 
  {
    return aCurrentCost > mThres ;
  }
  
private:
  
  FT mThres ;
};    

} //namespace internal

} //namespace CGAL

/// @endcond

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_MINIMUM_LENGTH_PREDICATE_H //
// EOF //
 
