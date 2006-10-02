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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_SET_COST_CACHE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_SET_COST_CACHE_H

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Cost_cache.h>

CGAL_BEGIN_NAMESPACE

namespace Surface_mesh_simplification
{

template<class ECM_, class GetCost_>    
class Set_cost_cache
{
public:

  typedef ECM_     ECM ;
  typedef GetCost_ GetCost ;
  
  typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
  
  typedef Cost_cache<ECM> Cache ;
  
  typedef typename Cache::Optional_cost_type Optional_cost_type ;
  
  typedef typename GetCost::Params CostParams ;
  typedef void                     PlacementParams ;

public :

  Set_cost_cache ( GetCost const& aGetCost ) : get_cost(aGetCost) {}

  void operator() ( Cache&                  rCache
                  , edge_descriptor const&  aEdge
                  , ECM&                    aSurface
                  , CostParams const*       aCostParams 
                  , PlacementParams const* 
                  ) const 
  {
    CGAL_assertion(aCostParams);
    CGAL_assertion( handle_assigned(aEdge) );

    Optional_cost_type lCost = get_cost(aEdge,aSurface,rCache,aCostParams);
    
    rCache = Cache(lCost);
    
  }                         
  
  GetCost get_cost ;
};    

} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_SET_COST_CACHE_H
// EOF //
 
