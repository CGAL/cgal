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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_SET_FULL_CACHE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_SET_FULL_CACHE_H

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Cost_and_placement_cache.h>


CGAL_BEGIN_NAMESPACE

namespace Surface_mesh_simplification
{

template<class ECM_, class SetPlan_>    
class Set_full_cache
{
public:

  typedef ECM_     ECM ;
  typedef SetPlan_ SetPlan ;
  
  typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
  
  typedef typename SetPlan::Params SetPlanParams ;
  typedef typename SetPlan::Plan   Plan ;
  
  typedef Plan Cache ;
  
public :

  Set_full_cache ( SetPlan const& aSetPlan ) : set_plan(aSetPlan) {}
    
  void operator() ( Cache&                 rCache
                  , edge_descriptor const& aEdge
                  , ECM&                   aSurface
                  , SetPlanParams const*   aSetPlanParams 
                  ) const 
  {
    CGAL_assertion( handle_assigned(aEdge) );

    set_plan(rCache,aEdge,aSurface,aSetPlanParams);    
  }                         
  
private :
  
  SetPlan set_plan ;
};    

} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_COST_AND_PLACEMET_CACHE_H
// EOF //
 
