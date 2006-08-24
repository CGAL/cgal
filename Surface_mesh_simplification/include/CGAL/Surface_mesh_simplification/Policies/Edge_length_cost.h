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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_LENGHT_COST_MAP_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_LENGHT_COST_MAP_H 1

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{


//*******************************************************************************************************************
//                               -= cost property maps =-
//
// Computes the cost of collapsing a vertex-pair.
// The map can return an empty optional if the cost is too-high or uncomputable (due to overflow for example).
//
//*******************************************************************************************************************

//
// Edge-length cost: the square distance between the collapsing vertices.
//
template<class Collapse_data_>
class Edge_length_cost
{
public:
    
  typedef Collapse_data_ Collapse_data ;
  
  typedef typename Collapse_data::TSM               TSM ;
  typedef typename Collapse_data::vertex_descriptor vertex_descriptor ;
  typedef typename Collapse_data::edge_descriptor   edge_descriptor ;
  typedef typename Collapse_data::Point_3           Point_3 ;
  typedef typename Collapse_data::FT                FT ;
    
public:

  typedef optional<FT> result_type;
    
  result_type operator()( edge_descriptor const& edge, TSM const& aSurface, Collapse_data const& data ) const
  {
    if ( !data.is_edge_fixed() )
    {
      Point_3 const& p = get_point(data.p(),data.surface());
      Point_3 const& q = get_point(data.q(),data.surface());
      
      return result_type(squared_distance(p,q));
    }
    else return result_type();
  }
  
  result_type compute_cost( edge_descriptor const& edge, TSM const& aSurface, Collapse_data const& data ) const
  {
    Point_3 const& p = get_point(data.p(),data.surface());
    Point_3 const& q = get_point(data.q(),data.surface());
      
    return result_type(squared_distance(p,q));
  }
  
  Point_3 const& get_point ( vertex_descriptor const& v, TSM const& aSurface ) const
  {
    vertex_point_t vertex_point ;
    return get(vertex_point,aSurface,v) ;
  }
  
  

};


} } // namespace Triangulated_surface_mesh::Simplification


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_LENGHT_COST_MAP_H //
// EOF //
 
