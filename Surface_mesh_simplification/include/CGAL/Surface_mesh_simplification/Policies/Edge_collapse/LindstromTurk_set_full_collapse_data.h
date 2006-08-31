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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_SET_FULL_COLLAPSE_DATA_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_SET_FULL_COLLAPSE_DATA_H 1

#include <CGAL/Surface_mesh_simplification/Detail/TSMS_common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Full_collapse_data.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_params.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core.h>


CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse
{

template<class TSM_>    
class Set_full_collapse_data_LindstromTurk
{
public:

  typedef TSM_ TSM ;

  typedef LindstromTurk_params Params ;
  
  typedef typename boost::graph_traits<TSM>::vertex_descriptor vertex_descriptor ;
  typedef typename boost::graph_traits<TSM>::edge_descriptor   edge_descriptor ;
  
  typedef typename Surface_geometric_traits<TSM>::FT      FT ;
  typedef typename Surface_geometric_traits<TSM>::Point_3 Point_3 ;
  
  typedef Full_collapse_data<TSM> Collapse_data ;
  
  typedef optional<FT>      Optional_cost_type ;
  typedef optional<Point_3> Optional_placement_type ;
  
  typedef LindstromTurk_params CostParams      ;
  typedef LindstromTurk_params PlacementParams ;

public :

  Set_full_collapse_data_LindstromTurk() {}
  
  void operator() ( Collapse_data&         rData
                  , edge_descriptor const& aEdge
                  , TSM&                   aSurface
                  , CostParams const*      aCostParams 
                  , PlacementParams const* aPlacementParams 
                  ) const 
  {
    CGAL_assertion(aCostParams);
    CGAL_assertion(aPlacementParams);
    CGAL_assertion(aCostParams == aPlacementParams);
    CGAL_assertion( handle_assigned(aEdge) );
    
    LindstromTurkCore<TSM> core(*aCostParams,aEdge,aSurface,true);

    Optional_cost_type lCost ;
    Optional_placement_type lPlacement ;
    tie(lCost,lPlacement) = core.compute();
    rData = Collapse_data(lCost,lPlacement);
  }                         
  
};    

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_SET_FULL_COLLAPSE_DATA_H //
// EOF //
 
