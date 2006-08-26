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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse
{

//
// Edge-length cost: the squared length of the collapsing edge
//
template<class TSM_>
class Edge_length_cost
{
public:
    
  typedef TSM_ TSM ;
  
  typedef typename boost::graph_traits<TSM>::edge_descriptor   edge_descriptor ;
  typedef typename boost::graph_traits<TSM>::vertex_descriptor vertex_descriptor ;
  
  typedef typename Surface_geometric_traits<TSM>::Point_3 Point_3 ;
  typedef typename Surface_geometric_traits<TSM>::FT      FT ;
  
  typedef optional<FT> result_type ;
    
public:

  mutable int c ; Edge_length_cost() : c(0) {}
  
  
  template<class Collapse_data, class Params>
  result_type operator()( edge_descriptor const& aEdge
                        , TSM&                   aSurface
                        , Collapse_data const&   aData
                        , Params const*          aParams
                        ) const
  {
    vertex_descriptor vs,vt ; tie(vs,vt) = this->get_vertices(aEdge,aSurface);
    
    Point_3 const& ps = this->get_point(vs,aSurface);
    Point_3 const& pt = this->get_point(vt,aSurface);
      
    ++ c ;
    
    return result_type(squared_distance(ps,pt));
  }
  
private:

  Point_3 const& get_point ( vertex_descriptor const& v, TSM& aSurface ) const
  {
    vertex_point_t vertex_point ;
    return get(vertex_point,aSurface,v) ;
  }
  
  tuple<vertex_descriptor,vertex_descriptor> get_vertices ( edge_descriptor const& aEdge, TSM& aSurface ) const
  {
    vertex_descriptor p = source(aEdge,aSurface);
    vertex_descriptor q = target(aEdge,aSurface);
    return make_tuple(p,q);
  }

};


} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H
// EOF //
 
