// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
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
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_CONNECTED_COMPONENTS_H
#define CGAL_VORONOI_DIAGRAM_2_CONNECTED_COMPONENTS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//========================================================================
//========================================================================

template<class VD_t>
class Connected_components
{
 private:
  typedef VD_t                                            VD;
  typedef typename VD::Halfedge_iterator                  Halfedge_iterator;
  typedef typename VD::Halfedge_handle                    Halfedge_handle;
  typedef typename VD::Halfedge                           Halfedge;
  typedef typename VD::Halfedge_around_vertex_circulator  HAVC;

  typedef HAVC Halfedge_around_vertex_circulator;

 public:
  typedef VD  Voronoi_diagram_2;

  typedef typename Voronoi_diagram_2::size_type          size_type;
  typedef size_type                                      result_type;
  typedef Voronoi_diagram_2                              argument_type;

 private:
  struct Halfedge_handle_less {
    bool operator()(const Halfedge_handle& e1,
		    const Halfedge_handle& e2) const {
      typename Halfedge::Delaunay_edge de1 = e1->dual();
      typename Halfedge::Delaunay_edge de2 = e2->dual();

      if ( de1.first != de2.first ) { return de1.first < de2.first; }
      return de1.second < de2.second;
    }
  };

  typedef std::map<Halfedge_handle,bool,Halfedge_handle_less>
  Halfedge_handle_map;

  void mark(const Halfedge_handle& e, Halfedge_handle_map& e_map) const
  {
    e_map[e] = true;
    e_map[e->opposite()] = true;
  }

  bool is_unmarked(const Halfedge_handle& e,
		   const Halfedge_handle_map& e_map) const
  {
    return e_map.find(e) == e_map.end();
  }

  void dfs(const Voronoi_diagram_2& vd, const Halfedge_handle& e,
	   Halfedge_handle_map& e_map) const
  {
    CGAL_precondition( !vd.dual().is_infinite(e->dual()) );

    Halfedge_handle e_opp = e->opposite();
    mark(e, e_map);

    if ( e->has_source() ) {
      HAVC ec =	vd.incident_halfedges(e->source());
      HAVC ec_start = ec;

      do {
	if ( e != ec && e_opp != ec && is_unmarked(ec, e_map) ) {
	  dfs(vd, ec, e_map);
	}
	ec++;
      } while (ec != ec_start);
    }

    if ( e->has_target() ) {
      HAVC ec = vd.incident_halfedges(e->target());
      HAVC ec_start = ec;

      do {
	if ( e != ec && e_opp != ec && is_unmarked(ec, e_map) ) {
	  dfs(vd, ec, e_map);
	}
	ec++;
      } while (ec != ec_start);
    }
  }

 public:
  size_type operator()(const Voronoi_diagram_2& vd) const
  {
    Halfedge_handle_map e_map;

    size_type n_components = 0;
    for(Halfedge_iterator eit = vd.halfedges_begin();
	eit != vd.halfedges_end(); ++eit) {
      if ( is_unmarked(eit, e_map) ) {
	n_components++;
	dfs(vd, eit, e_map);
      }
    }
    return n_components;
  }
};

//========================================================================
//========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_CONNECTED_COMPONENTS
