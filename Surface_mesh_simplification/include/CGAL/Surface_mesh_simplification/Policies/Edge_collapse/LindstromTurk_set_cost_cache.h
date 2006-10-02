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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_set_collapse_data.h $
// $Id: LindstromTurk_set_collapse_data.h 33680 2006-08-24 15:22:37Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_SET_COST_CACHE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_SET_COST_CACHE_H 1

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Cost_cache.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_params.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core.h>


CGAL_BEGIN_NAMESPACE

namespace Surface_mesh_simplification
{

template<class ECM_>    
class LindstromTurk_set_cost_cache
{
public:

  typedef ECM_ ECM ;

  typedef LindstromTurk_params Params ;
  
  typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
  
  typedef Cost_cache<ECM> Cache ;
  
  typedef typename halfedge_graph_traits<ECM>::Point Point ;
  typedef typename Kernel_traits<Point>::Kernel      Kernel ;
  typedef typename Kernel::FT                        FT ;
  
  typedef optional<FT>    Optional_cost_type ;
  typedef optional<Point> Optional_placement_type ;

  typedef LindstromTurk_params CostParams ;
  typedef void                 PlacementParams ;
  
public :

  LindstromTurk_set_cost_cache() {}
    
  void operator() ( Cache&                 rCache
                  , edge_descriptor const& aEdge
                  , ECM&                   aSurface
                  , CostParams const*      aCostParams 
                  , PlacementParams const* 
                  ) const 
  {
    CGAL_assertion(aCostParams);
    CGAL_assertion( handle_assigned(aEdge) );
    
    LindstromTurkCore<ECM> core(*aCostParams,aEdge,aSurface,true);

    Optional_cost_type lCost ;
    Optional_placement_type lPlacement ;
    tie(lCost,lPlacement) = core.compute();
    rCache = Cache(lCost);
}                         
  
};    

} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_SET_COST_CACHE_H //
// EOF //
 
