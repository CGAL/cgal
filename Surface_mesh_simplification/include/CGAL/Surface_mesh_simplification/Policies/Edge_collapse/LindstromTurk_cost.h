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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_COST_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_COST_H 1

#include <CGAL/Surface_mesh_simplification/Detail/ECMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Surface_mesh_simplification  
{

template<class ECM_>
class LindstromTurk_cost
{
public:
    
  typedef ECM_ ECM ;
  
  typedef typename boost::graph_traits<ECM>::vertex_descriptor vertex_descriptor ;
  typedef typename boost::graph_traits<ECM>::edge_descriptor   edge_descriptor ;
  
  typedef typename halfedge_graph_traits<ECM>::Point Point_3 ;
  typedef typename Kernel_traits<Point_3>::Kernel     Kernel ;
  typedef typename Kernel::FT                         FT ;
  
  typedef LindstromTurk_params Params ;
    
  typedef optional<FT> result_type ;
  
public:

  LindstromTurk_cost() {}
     
  template<class Cache>
  result_type operator()( edge_descriptor const& aEdge
                        , ECM&                   aSurface
                        , Cache const&           aCache
                        , Params const*          aParams
                        ) const
  {
    CGAL_assertion(aParams);
    CGAL_assertion( handle_assigned(aEdge) );
    
    LindstromTurkCore<ECM> core(*aParams,aEdge,aSurface,true);

    optional<FT> lCost ;
    optional<Point_3> lPlacement ;
    tie(lCost,lPlacement) = core.compute();
    
    return lCost ;
  }
  
};

} // namespace Surface_mesh_simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_COST_H //
// EOF //
 
