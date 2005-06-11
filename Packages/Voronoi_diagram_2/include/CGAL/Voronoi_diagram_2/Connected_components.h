// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_CONNECTED_COMPONENTS_H
#define CGAL_CONNECTED_COMPONENTS_H 1

template<class VDA_t>
class Connected_components
{
 private:
  typedef VDA_t                                            VDA;
  typedef typename VDA::Halfedge_iterator                  Halfedge_iterator;
  typedef typename VDA::Halfedge_handle                    Halfedge_handle;
  typedef typename VDA::Vertex_handle                      Vertex_handle;
  typedef typename VDA::Halfedge                           Halfedge;
  typedef typename VDA::Vertex                             Vertex;
  typedef typename VDA::Halfedge_around_vertex_circulator  HAVC;

  typedef HAVC Halfedge_around_vertex_circulator;

 public:
  typedef VDA  Voronoi_diagram_adaptor_2;

  typedef typename Voronoi_diagram_adaptor_2::size_type  size_type;
  typedef size_type                                      result_type;
  typedef Voronoi_diagram_adaptor_2                      argument_type;

  struct Arity { enum { arity = 1 }; };

 private:
  struct Halfedge_less {
    bool operator()(const Halfedge& h1,	const Halfedge& h2) const {
      typename Voronoi_diagram_adaptor_2::Dual_edge e1 = h1.dual_edge();
      typename Voronoi_diagram_adaptor_2::Dual_edge e2 = h2.dual_edge();

      if ( e1.first != e2.first ) { return e1.first < e2.first; }
      return e1.second < e2.second;
    }
  };

  typedef std::map<Halfedge,bool,Halfedge_less>  Halfedge_map;

  void mark_edge(const Halfedge& e, Halfedge_map& e_map) const
  {
    e_map[e] = true;
    e_map[*e.opposite()] = true;
  }

  bool is_unmarked(const Halfedge& e, const Halfedge_map& e_map) const
  {
    return e_map.find(e) == e_map.end();
  }

  void dfs(const Voronoi_diagram_adaptor_2& vda, const Halfedge& e,
	      Halfedge_map& e_map) const
  {
    CGAL_precondition( !vda.dual().is_infinite(e.dual_edge()) );

    Halfedge e_opp = *e.opposite();
    mark_edge(e, e_map);

    if ( e.has_source() ) {
      HAVC ec =	vda.incident_halfedges(e.source());
      HAVC ec_start = ec;

      do {
	if ( *ec != e && *ec != e_opp && is_unmarked(*ec, e_map) ) {
	  dfs(vda, *ec, e_map);
	}
	ec++;
      } while (ec != ec_start);
    }

    if ( e.has_target() ) {
      HAVC ec = vda.incident_halfedges(e.target());
      HAVC ec_start = ec;

      do {
	if ( *ec != e && *ec != e_opp && is_unmarked(*ec, e_map) ) {
	  dfs(vda, *ec, e_map);
	}
	ec++;
      } while (ec != ec_start);
    }
  }

 public:
  size_type operator()(const Voronoi_diagram_adaptor_2& vda) const
  {
    Halfedge_map e_map;

    size_type n_components = 0;
    for(Halfedge_iterator eit = vda.halfedges_begin();
	eit != vda.halfedges_end(); ++eit) {
      if ( is_unmarked(*eit, e_map) ) {
	n_components++;
	dfs(vda, *eit, e_map);
      }
    }
    return n_components;
  }
};


#endif // CGAL_CONNECTED_COMPONENTS
