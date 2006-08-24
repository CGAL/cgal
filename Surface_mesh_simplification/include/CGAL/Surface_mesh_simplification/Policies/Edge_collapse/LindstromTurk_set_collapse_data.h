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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_SET_COLLAPSE_DATA_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_SET_COLLAPSE_DATA_H 1

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>
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
  
  typedef typename boost::graph_traits<TSM>::edge_descriptor edge_descriptor ;
  
  typedef Full_collapse_data<TSM> Collapse_data ;
  
  typedef typename Collapse_data::Optional_cost_type      Optional_cost_type ;
  typedef typename Collapse_data::Optional_placement_type Optional_placement_type ;

public :

  void operator() ( Collapse_data& rData, edge_descriptor const& aEdge, TSM& aSurface, Params const* aParams ) const 
  {
    CGAL_assertion(aParams);
    CGAL_assertion( handle_assigned(aEdge) );
    
    LindstromTurkCore<TSM> core(*aParams,aEdge,aSurface,true);

    Optional_cost_type lCost ;
    Optional_placement_type lPlacement ;
    tie(lCost,lPlacement) = core.compute();
    rData.set(lCost,lPlacement);
  }                         
  
};    

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_SET_COLLAPSE_DATA_H //
// EOF //
 
