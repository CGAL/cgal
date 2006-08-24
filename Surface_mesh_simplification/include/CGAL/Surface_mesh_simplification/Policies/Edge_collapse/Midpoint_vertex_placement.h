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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_MIDPOINT_VERTEX_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_MIDPOINT_VERTEX_PLACEMENT_H 1

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse
{

//
// Edge-length cost: the squared length of the collapsing edge
//
template<class Collapse_data_>
class Midpoint_placement : protected Placement_functor_base< Collapse_data_, Midpoint_placement<Collapse_data_> >
{
  typedef Cost_functor_base<Collapse_data_, Midpoint_placement<Collapse_data_> > Base ;
  
public:
    
  typedef Collapse_data_ Collapse_data ;
  
  typedef typename Collapse_data::TSM TSM ;
  
  typedef typename grah_traits<TSM>::vertex_descriptor vertex_descriptor ;
  typedef typename grah_traits<TSM>::edge_descriptor   edge_descriptor ;
  
  typedef typename Surface_geometric_traits<TSM>::Point_3 Point_3 ;
  
  typedef typename Base::result_type result_type ;
    
public:

  result_type compute_placement( edge_descriptor const& aEdge, TSM const& aSurface ) const
  {
    vertex_descriptor vs,vt ; tie(vs,vt) = this->get_vertices(aEdge,aSurface);
    
    Point_3 const& ps = this->get_point(vs,aSurface);
    Point_3 const& pt = this->get_point(vt,aSurface);
      
    return result_type(midpoint(ps,pt));
  }
};

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_MIDPOINT_VERTEX_PLACEMENT_H //
// EOF //
 
