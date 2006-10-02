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

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>

CGAL_BEGIN_NAMESPACE

namespace Surface_mesh_simplification
{

//
// Edge-length cost: the squared length of the collapsing edge
//
template<class ECM_>
class Edge_length_cost
{
  
public:
    
  typedef ECM_ ECM ;
  
private :
    
  typedef typename boost::graph_traits<ECM>::vertex_descriptor vertex_descriptor ;
  
  typedef typename halfedge_graph_traits<ECM>::Point Point ;
  
  typedef typename Kernel_traits<Point>::Kernel Kernel ;

public:

  typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
  
  typedef typename Kernel::FT FT ;
  
  typedef optional<FT> result_type ;
  
  typedef char Params ;

  
public:

  Edge_length_cost() {}
  
  template<class Cache>
  result_type operator()( edge_descriptor const& aEdge
                        , ECM&                   aSurface
                        , Cache const&           aCache
                        , Params const*          aParams
                        ) const
  {
    vertex_descriptor vs,vt ; tie(vs,vt) = this->get_vertices(aEdge,aSurface);
    
    Point_3 const& ps = get(vertex_point,aSurface,vs);
    Point_3 const& pt = get(vertex_point,aSurface,vt);
      
    return result_type(squared_distance(ps,pt));
  }
  
private:
  
  tuple<vertex_descriptor,vertex_descriptor> get_vertices ( edge_descriptor const& aEdge, ECM& aSurface ) const
  {
    vertex_descriptor p = source(aEdge,aSurface);
    vertex_descriptor q = target(aEdge,aSurface);
    return make_tuple(p,q);
  }

};


} // namespace Surface_mesh_simplification


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H
// EOF //
 
