// Copyright (c) 2017  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_NORMAL_CHANGE_PLACEMENT_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_NORMAL_CHANGE_PLACEMENT_H

#include <CGAL/license/Surface_mesh_simplification.h>
#include <boost/optional.hpp>

namespace CGAL {

namespace Surface_mesh_simplification  
{

template<class Placement>
class Bounded_normal_change_placement
{
public:
    
  typedef typename Placement::ECM ECM ;
  
public:
  
  Bounded_normal_change_placement(const Placement& placement = Placement() )
    : mPlacement(placement)
  {}
     
  template <typename Profile> 
  boost::optional<typename Profile::Point>
  operator()( Profile const& aProfile) const
  {
    boost::optional<typename Profile::Point> op = mPlacement(aProfile);
    if(op){
      // triangles returns the triangles of the star of the vertices of the edge to collapse
      // First the two trianges incident to the edge, then the other triangles
      // The second vertex of each triangle is the vertex that gets placed
       const typename Profile::Triangle_vector& triangles = aProfile.triangles();
       if(triangles.size()>2){
         typedef typename Profile::Point Point;
         typedef typename Profile::Kernel Traits;
         typedef typename Traits::Vector_3 Vector;
         typename Profile::VertexPointMap ppmap = aProfile.vertex_point_map();
         typename Profile::Triangle_vector::const_iterator it = triangles.begin();
         if(aProfile.left_face_exists()){
           ++it; 
         }
         if(aProfile.right_face_exists()){
           ++it;
         }
         while(it!= triangles.end()){
           const typename Profile::Triangle& t = *it;
           Point p = get(ppmap,t.v0);
           Point q = get(ppmap,t.v1);
           Point r = get(ppmap,t.v2);
           Point q2 = *op;
           
           Vector eqp = Traits().construct_vector_3_object()(q,p) ;
           Vector eqr = Traits().construct_vector_3_object()(q,r) ;
           Vector eq2p = Traits().construct_vector_3_object()(q2,p) ;
           Vector eq2r = Traits().construct_vector_3_object()(q2,r) ;
           
           Vector n1 = Traits().construct_cross_product_vector_3_object()(eqp,eqr);
           Vector n2 = Traits().construct_cross_product_vector_3_object()(eq2p,eq2r);
           if(! is_positive(Traits().compute_scalar_product_3_object()(n1, n2))){
             return boost::optional<typename Profile::Point>();
           }
           ++it;
         }
       }
    }
    return op;
  }
  
private:

  Placement  mPlacement ;

};


} // namespace Surface_mesh_simplification

} //namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_BOUNDED_NORMAL_CHANGE_PLACEMENT_H

 
