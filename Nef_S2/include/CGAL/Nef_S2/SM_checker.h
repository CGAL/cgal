// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// 
//
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_SM_CHECKER_H
#define CGAL_SM_CHECKER_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/basic.h>
#include <CGAL/Nef_S2/SM_const_decorator.h>

#define CGAL_USING(t) typedef typename Base::t t
namespace CGAL {

/*{\Moptions outfile=SM_checker.man }*/
/*{\Manpage {SM_checker}{PMCDEC,GEOM}{Plane map checking}{}}*/

/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
decorator to check the structure of a plane map. It is generic with
respect to two template concepts.  |PMCDEC| has to be a decorator
model of our |SM_const_decorator| concept. |GEOM| has to be a model of
our geometry kernel concept.}*/

/*{\Mgeneralization SM_const_decorator}*/

template <typename PMCDEC, typename GEOM> 
class SM_checker : public PMCDEC
{ typedef PMCDEC Base;
  const GEOM& K;
public:
/*{\Mtypes 3}*/
typedef PMCDEC  SM_const_decorator;
/*{\Mtypemember equals |PMCDEC|. Add link to PMCDEC concept.}*/
typedef typename PMCDEC::Plane_map Plane_map;
/*{\Mtypemember equals |PMCDEC::Plane_map|, the underlying plane map type.}*/
typedef GEOM Geometry;
/*{\Mtypemember equals |GEOM|. Add link to GEOM concept.\\
\precond |Geometry::Point_2| equals |Plane_map::Point|. }*/

typedef typename GEOM::Point_2     Point;
typedef typename GEOM::Direction_2 Direction;

CGAL_USING(Vertex_const_handle);
CGAL_USING(SEdge_handle);
CGAL_USING(Vertex_const_iterator);
CGAL_USING(Halfedge_const_iterator);
CGAL_USING(Halfedge_around_vertex_const_circulator);
CGAL_USING(Halfedge_around_face_const_circulator);
/*{\Mtext Iterators, handles, and circulators are inherited from 
|SM_const_decorator|.}*/

/*{\Mcreation 3}*/
SM_checker(Plane_map& P, const Geometry& k = Geometry()) : 
  Base(P), K(k) {}
/*{\Mcreate constructs a plane map checker working on |P| with
geometric predicates used from |k|.}*/

SM_checker(const Base& D, const Geometry& k = Geometry()) : 
  Base(D), K(k) {}


/*{\Moperations 2 }*/
Direction direction(Halfedge_const_handle e) const
{ return K.construct_direction(
    point(source(e)),point(target(e))); }

bool is_forward(Halfedge_const_handle e) const
{ return K.compare_xy(point(source(e)),point(target(e)))<0; }

void check_order_preserving_embedding(Vertex_const_handle v) const;
/*{\Mop checks if the embedding of the targets of the edges in
the adjacency list |A(v)| is counter-clockwise order-preserving with 
respect to the order of the edges in |A(v)|.}*/

void check_order_preserving_embedding() const;
/*{\Mop checks if the embedding of all vertices of |P| is 
counter-clockwise order-preserving with respect to the adjacency
list ordering of all vertices.}*/

void check_forward_prefix_condition(Vertex_const_handle v) const;
/*{\Mop checks the forward prefix property of the adjacency
list of |v|.}*/

Halfedge_const_iterator
check_boundary_is_clockwise_weakly_polygon() const;
/*{\Mop checks if the outer face cycle of |P| is a clockwise weakly polygon and
returns a halfedge on the boundary. \precond |P| is a connected graph.
}*/

void check_is_triangulation() const;
/*{\Mop checks if |P| is a triangulation.}*/

}; // SM_checker<PMCDEC,GEOM>


template <typename PMCDEC, typename GEOM>
void SM_checker<PMCDEC,GEOM>::
check_order_preserving_embedding(Vertex_const_handle v) const
{
  std::ostrstream error_status;
  CGAL::set_pretty_mode ( error_status );
  Halfedge_const_handle ef = first_out_edge(v) ,e=ef,en,enn;
  error_status << "check_order_preserving_embedding\n";
  error_status << "vertex " << PV(v) << endl;
  if ( e != Halfedge_const_handle() ) {
    while ( true ) {
      en = cyclic_adj_succ(e);
      enn = cyclic_adj_succ(en);
      if (en == ef) break;
      error_status << "  -> " << point(target(e));
      error_status << " " << point(target(en)) << " ";
      error_status << " " << point(target(enn)) << endl;
      if ( !K.strictly_ordered_ccw(direction(e),direction(en),
                                   direction(enn)) ||
           !K.strictly_ordered_ccw(direction(e),direction(en),
                                   direction(ef)) ) {
        error_status << "ccw order violate!" << endl << '\0';
        CGAL_error_msg(error_status.str());
      }
      e = en;
    }
  }
  error_status.freeze(0);  
}

template <typename PMCDEC, typename GEOM>
void SM_checker<PMCDEC,GEOM>::
check_forward_prefix_condition(Vertex_const_handle v) const
{ Halfedge_const_handle ef = first_out_edge(v);
  if ( ef == Halfedge_const_handle() ) return;
  Halfedge_const_handle el = cyclic_adj_pred(ef);
  bool is_left_turn = K.left_turn(point(v),
                                point(target(ef)),
                                point(target(el)));
  bool el_forward = is_forward(el);
  bool ef_forward = is_forward(ef);
  CGAL_assertion_msg( (ef == el ||
                       ef_forward && !el_forward ||
                       ef_forward &&  el_forward && is_left_turn ||
                       !ef_forward && !el_forward && is_left_turn) ,
  "check_forward_prefix_condition: first backward, last forward\n");
}

template <typename PMCDEC, typename GEOM>
void SM_checker<PMCDEC,GEOM>::
check_order_preserving_embedding() const
{
  Vertex_const_iterator vit;
  for (vit = vertices_begin(); vit != vertices_end(); ++vit) {
    check_order_preserving_embedding(vit);
    check_forward_prefix_condition(vit);
  }
}


template <typename PMCDEC, typename GEOM>
typename SM_checker<PMCDEC,GEOM>::Halfedge_const_iterator
SM_checker<PMCDEC,GEOM>::
check_boundary_is_clockwise_weakly_polygon() const
{
  Vertex_const_iterator vit, v_min;
  for (vit = v_min = vertices_begin() ; vit != vertices_end(); ++vit) 
    if ( K.compare_xy(point(vit), point(v_min))<0 ) v_min = vit;
  CGAL_assertion_msg(!is_isolated(v_min),"Minimal vertex not connected.");
  Sphere_point p_min = point(v_min);

  // determine boundary edge incident to v_min: 
  Halfedge_const_handle e_boundary_at_v_min = first_out_edge(v_min);
  // all out edges are forward oriented due to minimality
  Halfedge_around_vertex_const_circulator 
    hvit(e_boundary_at_v_min), hend(hvit);
  do {
    --hvit;
    Sphere_point p1 = point(target(e_boundary_at_v_min));
    Sphere_point p2 = point(target(hvit));
    if ( K.orientation(p_min,p1,p2) > 0 ) { // left_turn
      e_boundary_at_v_min = hvit;
      break;
    }
  } while (hvit != hend);
  // now e_boundary_at_v_min is highest starting edge in bundle!!

  int winding_around_globally=0;
  Halfedge_around_face_const_circulator
    hfit(e_boundary_at_v_min),hstart(hfit);
  Halfedge_const_handle e_prev = next(e_boundary_at_v_min);
  /* we run counterclockwise around the outer face cycle and allow only
     position where the direction vector of an edge gets smaller again */
  Direction d_prev = direction(e_prev);
  CGAL_For_all_backwards(hstart,hfit) {
    Direction d_curr = direction(hfit);
    if ( d_curr < d_prev ) ++winding_around_globally;
    d_prev = d_curr;
  }
  CGAL_assertion(winding_around_globally == 1);
  return e_boundary_at_v_min;
}

template <typename PMCDEC, typename GEOM>
void SM_checker<PMCDEC,GEOM>::
check_is_triangulation() const
{
  Halfedge_const_iterator eb;
  check_integrity_and_topological_planarity(false);
  CGAL_assertion(number_of_connected_components() == 1);
  check_order_preserving_embedding();
  eb = check_boundary_is_clockwise_weakly_polygon();

  CGAL::Hash_map< Halfedge_const_iterator, bool> on_boundary(false);
  Halfedge_around_face_const_circulator hit(eb), hend(hit);
  std::ostrstream error_status;
  CGAL::set_pretty_mode ( error_status );
  error_status << "check_is_triangulation\n";
  error_status << "on boundary:\n";
  CGAL_For_all(hit,hend) {
    error_status << "  " << PE(hit) << endl;
    on_boundary[hit]=true;
  }
  Halfedge_const_iterator eit;
  for( eit = halfedges_begin(); eit != halfedges_end(); ++eit) {
    if (on_boundary[eit]) continue;
    hit = hend = eit; 
    int edges_in_face_cycle=0;
    CGAL_For_all(hit,hend) {
      error_status << PE(hit);
      ++edges_in_face_cycle;
    }
    CGAL_assertion_msg(edges_in_face_cycle==3,error_status.str());
  }
  error_status.freeze(0);
}


/*\subsection{Plane map input and output}
The input and output is mainly triggered by an IO Decorator which
has the control over the IO format and does some basic parsing when
reading input.*/

} //namespace CGAL

#undef CGAL_USING
#endif // CGAL_SM_CHECKER_H
