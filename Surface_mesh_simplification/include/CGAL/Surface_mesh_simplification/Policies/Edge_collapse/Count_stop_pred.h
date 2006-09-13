// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_COUNT_STOP_PRED_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_COUNT_STOP_PRED_H 1

#include <CGAL/Surface_mesh_simplification/Detail/TSMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse
{

//*******************************************************************************************************************
//                                -= stopping condition predicate =-
//
// Determines whether the simplification has finished.
// The arguments are (current_cost,vertex,vertex,is_edge,initial_pair_count,current_pair_count,surface) and the result is bool
//
//*******************************************************************************************************************

// 
// Stops when the number of edges left falls below a given number.
//
template<class TSM_>    
class Count_stop_condition
{
public:

  typedef TSM_ TSM ;
  
  typedef typename boost::graph_traits<TSM>::edge_descriptor edge_descriptor ;
  typedef typename boost::graph_traits<TSM>::edges_size_type size_type ;
  
  typedef typename Geometric_graph_traits<TSM>::Point Point_3 ;
  typedef typename Kernel_traits<Point_3>::Kernel     Kernel ;
  typedef typename Kernel::FT                         FT ;

public :
  
  Count_stop_condition( size_type aThres ) : mThres(aThres) {}
  
  bool operator()( FT const&              // aCurrentCost
                 , edge_descriptor const& //aEdge
                 , size_type              aInitialCount
                 , size_type              aCurrentCount
                 ) const 
  {
    return aCurrentCount < mThres ;
  }
  
private:
  
  size_type mThres ;
};    

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_COUNT_STOP_PRED_H //
// EOF //
 
