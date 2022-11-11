// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
