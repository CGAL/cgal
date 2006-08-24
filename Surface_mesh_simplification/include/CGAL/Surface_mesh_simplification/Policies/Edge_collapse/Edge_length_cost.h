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
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Functor_base.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse
{

//
// Edge-length cost: the squared length of the collapsing edge
//
template<class Collapse_data_>
class Edge_length_cost : public Cost_functor_base< Collapse_data_,void>
{
  typedef Cost_functor_base<Collapse_data_,void> Base ;
  
public:
    
  typedef Collapse_data_ Collapse_data ;
  
  typedef typename Collapse_data::TSM TSM ;
  
  typedef typename boost::graph_traits<TSM>::vertex_descriptor vertex_descriptor ;
  typedef typename boost::graph_traits<TSM>::edge_descriptor   edge_descriptor ;
  
  typedef typename Surface_geometric_traits<TSM>::Point_3 Point_3 ;
  
  typedef typename Base::result_type result_type ;
  typedef typename Base::Params      Params ;
    
public:

  virtual result_type compute_cost( edge_descriptor const& aEdge, TSM& aSurface, Params const* ) const
  {
    vertex_descriptor vs,vt ; tie(vs,vt) = this->get_vertices(aEdge,aSurface);
    
    Point_3 const& ps = this->get_point(vs,aSurface);
    Point_3 const& pt = this->get_point(vt,aSurface);
      
    return result_type(squared_distance(ps,pt));
  }
};


} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGHT_COST_H
// EOF //
 
