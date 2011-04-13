// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/PM_checker.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Geometric properties checker 
// ============================================================================

#ifndef CGAL_PM_CHECKER_H
#define CGAL_PM_CHECKER_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/PM_const_decorator.h>

#define USING(t) typedef typename Base::t t
CGAL_BEGIN_NAMESPACE

/*{\Moptions outfile=PM_checker.man }*/
/*{\Manpage {PM_checker}{PMCDEC,GEOM}{Plane map checking}{}}*/

/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
decorator to check the structure of a plane map. It is generic with
respect to two template concepts.  |PMCDEC| has to be a decorator
model of our |PM_const_decorator| concept. |GEOM| has to be a model of
our geometry kernel concept.}*/

/*{\Mgeneralization PM_const_decorator}*/

template <typename PMCDEC, typename GEOM> 
class PM_checker : public PMCDEC
{ typedef PMCDEC Base;
  const GEOM& K;
public:
/*{\Mtypes 3}*/
typedef PMCDEC  PM_const_decorator;
/*{\Mtypemember equals |PMCDEC|.}*/
typedef typename PMCDEC::Plane_map Plane_map;
/*{\Mtypemember equals |PMCDEC::Plane_map|, the underlying plane map type.}*/
typedef GEOM Geometry;
/*{\Mtypemember equals |GEOM|. Add link to GEOM concept.\\
\precond |Geometry::Point_2| equals |Plane_map::Point|. }*/

typedef typename GEOM::Point_2     Point;
typedef typename GEOM::Direction_2 Direction;

USING(Vertex_const_handle);
USING(Halfedge_const_handle);
USING(Vertex_const_iterator);
USING(Halfedge_const_iterator);
USING(Halfedge_around_vertex_const_circulator);
USING(Halfedge_around_face_const_circulator);
/*{\Mtext Iterators, handles, and circulators are inherited from 
|PM_const_decorator|.}*/

/*{\Mcreation 3}*/
PM_checker(Plane_map& P, const Geometry& k = Geometry()) : 
  Base(P), K(k) {}
/*{\Mcreate constructs a plane map checker working on |P| with
geometric predicates used from |k|.}*/

PM_checker(const Base& D, const Geometry& k = Geometry()) : 
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
/*{\Mop checks the forward-prefix property of the adjacency list of |v|.}*/

Halfedge_const_iterator
check_boundary_is_clockwise_weakly_polygon() const;
/*{\Mop checks if the outer face cycle of |P| is a clockwise weakly polygon 
and returns a halfedge on the boundary. \precond |P| is a connected graph.
}*/

void check_is_triangulation() const;
/*{\Mop checks if |P| is a triangulation.}*/

}; // PM_checker<PMCDEC,GEOM>


template <typename PMCDEC, typename GEOM>
void PM_checker<PMCDEC,GEOM>::
check_order_preserving_embedding(Vertex_const_handle v) const
{
  if ( is_isolated(v) ) return;
  std::ostrstream error_status;
  CGAL::set_pretty_mode ( error_status );
  Halfedge_const_handle ef = first_out_edge(v) ,e=ef,en,enn;
  error_status << "check_order_preserving_embedding\n";
  error_status << "vertex " << PV(v) << std::endl;
  error_status << "ef " << PE(ef) << std::endl;
  while ( true ) {
    en = cyclic_adj_succ(e);
    enn = cyclic_adj_succ(en);
    if (en == ef) break;
    error_status << "  -> " << point(target(e))
                 << " " << point(target(en))  
                 << " " << point(target(enn)) << std::endl;
    bool ccw1 = K.strictly_ordered_ccw(direction(e),direction(en),
				       direction(enn));
    bool ccw2 = K.strictly_ordered_ccw(direction(e),direction(en),
				       direction(ef));
    if ( !(ccw1 && ccw2) ) {
      error_status << "ccw order violate!" << std::endl << '\0';
      CGAL_assertion_msg(0,error_status.str());
    }
    e = en;
  }
  error_status.freeze(0);  
}

template <typename PMCDEC, typename GEOM>
void PM_checker<PMCDEC,GEOM>::
check_forward_prefix_condition(Vertex_const_handle v) const
{ Halfedge_const_handle ef = first_out_edge(v);
  if ( ef == Halfedge_const_handle() ) return;
  Halfedge_const_handle el = cyclic_adj_pred(ef);
  bool is_leftturn = K.leftturn(point(v),
                                point(target(ef)),
                                point(target(el)));
  bool el_forward = is_forward(el);
  bool ef_forward = is_forward(ef);
  bool ef_el_eq = (ef==el);
  std::ostrstream error_status;
  error_status << "check_forward_prefix_condition: ";
  error_status << PV(v) << "\n";
  error_status << PE(ef) << "\n" << PE(el) << "\n";
  error_status << " ef == el = " << ef_el_eq;
  error_status << " ef_forward = " << ef_forward;
  error_status << " el_forward = " << el_forward;
  error_status << " is_leftturn = " << is_leftturn << '\0';
  CGAL_assertion_msg( (ef == el ||
                       ef_forward && !el_forward ||
                       ef_forward &&  el_forward && is_leftturn ||
                       !ef_forward && !el_forward && is_leftturn) ,
                       error_status.str());
  error_status.freeze(0);  
}

/* We check the geometric integrity of the structure. We check
   + that all adjacent nodes are differently embedded
   + that all node lists are correctly embedded counterclockwise
     with winding number one.
   + that the convex hull of the structure has winding number one.
*/

template <typename PMCDEC, typename GEOM>
void PM_checker<PMCDEC,GEOM>::
check_order_preserving_embedding() const
{
  Vertex_const_iterator vit;
  for (vit = vertices_begin(); vit != vertices_end(); ++vit) {
    check_order_preserving_embedding(vit);
    check_forward_prefix_condition(vit);
  }
}


template <typename PMCDEC, typename GEOM>
typename PM_checker<PMCDEC,GEOM>::Halfedge_const_iterator
PM_checker<PMCDEC,GEOM>::
check_boundary_is_clockwise_weakly_polygon() const
{
  Vertex_const_iterator vit, v_min;
  for (vit = v_min = vertices_begin() ; vit != vertices_end(); ++vit) 
    if ( K.compare_xy(point(vit), point(v_min))<0 ) v_min = vit;
  CGAL_assertion_msg(!is_isolated(v_min),"Minimal vertex not connected.");
  Point p_min = point(v_min);
  // determine boundary edge incident to v_min: 
  Halfedge_const_handle e_boundary_at_v_min = first_out_edge(v_min);
  // all out edges are forward oriented due to minimality
  Halfedge_around_vertex_const_circulator 
    hvit(e_boundary_at_v_min), hend(hvit);
  do {
    --hvit;
    Point p1 = point(target(e_boundary_at_v_min));
    Point p2 = point(target(hvit));
    if ( K.orientation(p_min,p1,p2) > 0 ) { // leftturn
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
void PM_checker<PMCDEC,GEOM>::
check_is_triangulation() const
{
  Halfedge_const_iterator eb;
  CGAL_assertion(number_of_connected_components() == 1);
  CGAL_assertion_msg(number_of_edges()!=number_of_vertices()-1,
    " checker checks only full dimensional complexes.");
  check_integrity_and_topological_planarity(false);
  check_order_preserving_embedding();
  eb = check_boundary_is_clockwise_weakly_polygon();

  CGAL::Unique_hash_map< Halfedge_const_iterator, bool> on_boundary(false);
  Halfedge_around_face_const_circulator hit(eb), hend(hit);
  std::ostrstream error_status;
  CGAL::set_pretty_mode ( error_status );
  error_status << "check_is_triangulation\n";
  error_status << "on boundary:\n";
  CGAL_For_all(hit,hend) {
    error_status << "  " << PE(hit) << std::endl;
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
    CGAL_assertion_msg(
      K.leftturn(point(source(hit)),point(target(hit)),
                 point(target(next(hit)))), error_status.str());
  }
  error_status.freeze(0);
}



CGAL_END_NAMESPACE

#undef USING
#endif // CGAL_PM_CHECKER_H

