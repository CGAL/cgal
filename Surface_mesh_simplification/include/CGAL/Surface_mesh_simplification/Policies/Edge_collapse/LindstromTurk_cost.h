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

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse  
{

template<class TSM_>
class LindstromTurk_cost
{
public:
    
  typedef TSM_ TSM ;
  
  typedef typename boost::graph_traits<TSM>::vertex_descriptor vertex_descriptor ;
  typedef typename boost::graph_traits<TSM>::edge_descriptor   edge_descriptor ;
  
  typedef typename Surface_geometric_traits<TSM>::FT      FT ;
  typedef typename Surface_geometric_traits<TSM>::Point_3 Point_3 ;
  
  typedef LindstromTurk_params Params ;
    
  typedef optional<FT> result_type ;
  
public:

  mutable int c ;

  LindstromTurk_cost() : c(0) {}
     
  template<class Collapse_data>
  result_type operator()( edge_descriptor const& aEdge
                        , TSM&                   aSurface
                        , Collapse_data const&   aData
                        , Params const*          aParams
                        ) const
  {
    CGAL_assertion(aParams);
    CGAL_assertion( handle_assigned(aEdge) );
    
    LindstromTurkCore<TSM> core(*aParams,aEdge,aSurface,true);

    optional<FT> lCost ;
    optional<Point_3> lPlacement ;
    tie(lCost,lPlacement) = core.compute();
    
    ++ c ;
    return lCost ;
  }
  
};

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_LINDSTROMTURK_COST_H //
// EOF //
 
