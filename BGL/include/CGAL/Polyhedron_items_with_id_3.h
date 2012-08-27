// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_POLYHEDRON_ITEMS_WITH_ID_3_H
#define CGAL_POLYHEDRON_ITEMS_WITH_ID_3_H 1

#include <CGAL/HalfedgeDS_vertex_max_base_with_id.h>
#include <CGAL/HalfedgeDS_halfedge_max_base_with_id.h>
#include <CGAL/HalfedgeDS_face_max_base_with_id.h>

namespace CGAL {

class Polyhedron_items_with_id_3 {
public:
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3 Point;
        typedef HalfedgeDS_vertex_max_base_with_id< Refs, Point, std::size_t> Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper {
        typedef HalfedgeDS_halfedge_max_base_with_id<Refs, std::size_t> Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef HalfedgeDS_face_max_base_with_id< Refs, Tag_false, std::size_t>  Face;
    };
};

template<class HalfedgeDS_with_id>
void set_halfedgeds_items_id ( HalfedgeDS_with_id& hds )
{
  std::size_t vertex_id   = 0 ;
  std::size_t halfedge_id = 0 ;
  std::size_t face_id     = 0 ;
  
  for ( typename HalfedgeDS_with_id::Vertex_iterator vit = hds.vertices_begin(), evit = hds.vertices_end()
      ; vit != evit
      ; ++  vit
      )
    vit->id() = vertex_id ++ ;    
    
  for ( typename HalfedgeDS_with_id::Halfedge_iterator hit = hds.halfedges_begin(), ehit = hds.halfedges_end()
      ; hit != ehit
      ; ++  hit
      )
    hit->id() = halfedge_id ++ ;    
    
  for ( typename HalfedgeDS_with_id::Face_iterator fit = hds.facets_begin(), efit = hds.facets_end()
      ; fit != efit
      ; ++  fit
      )
    fit->id() = face_id ++ ;    
}

} //namespace CGAL

#endif // CGAL_POLYHEDRON_ITEMS_WITH_ID_3_H //
// EOF //
