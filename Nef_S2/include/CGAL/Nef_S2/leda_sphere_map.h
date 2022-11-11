// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_LEDA_SPHERE_MAP_H
#define CGAL_LEDA_SPHERE_MAP_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/generic_sweep.h>
#include <CGAL/Nef_2/Segment_overlay_traits.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <LEDA/graph/graph.h>
#include <LEDA/graph/graph_misc.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 211
#include <CGAL/Nef_2/debug.h>

template <typename R, typename ITERATOR>
class leda_graph_decorator {
public:
  typedef leda_node      Vertex_handle;
  typedef leda_edge      Halfedge_handle;
  typedef CGAL::Sphere_point<R>   Point_2;
  typedef CGAL::Sphere_segment<R> Segment_2;
  typedef GRAPH< Point_2, Segment_2 >     Graph;
  typedef leda_node_map<Halfedge_handle>  Below_map;

  Graph&      G;
  Below_map&  M;

leda_graph_decorator(Graph& Gi, Below_map& Mi) : G(Gi), M(Mi) {}

Vertex_handle new_vertex(const Point_2& p)
{ return G.new_node(p); }

void link_as_target_and_append(Vertex_handle v, Halfedge_handle e)
{ Halfedge_handle erev = G.reversal(e);
  G.move_edge(e,G.cyclic_adj_pred(e,G.source(e)),v);
  G.move_edge(erev,v,G.target(erev));
}

Halfedge_handle new_halfedge_pair_at_source(Vertex_handle v)
{ Halfedge_handle e_res,e_rev, e_first = G.first_adj_edge(v);
  if ( e_first == nil ) {
    e_res = G.new_edge(v,v);
    e_rev = G.new_edge(v,v);
  } else {
    e_res = G.new_edge(e_first,v,LEDA::before);
    e_rev = G.new_edge(e_first,v,LEDA::before);
  }
  G.set_reversal(e_res,e_rev);
  return e_res;
}

void supporting_segment(Halfedge_handle e, ITERATOR it)
{ G[e] = *it; }

void halfedge_below(Vertex_handle v, Halfedge_handle e)
{ M[v] = e; }

void trivial_segment(Vertex_handle v, ITERATOR it) {}
void starting_segment(Vertex_handle v, ITERATOR it) {}
void passing_segment(Vertex_handle v, ITERATOR it) {}
void ending_segment(Vertex_handle v, ITERATOR it) {}


}; // leda_graph_decorator


template <typename R>
class leda_sphere_map_overlayer {

  typedef std::pair<leda_edge,leda_edge> edge_pair;
  typedef CGAL::Sphere_point<R>   SPoint_2;
  typedef CGAL::Sphere_segment<R> SSegment_2;
  typedef CGAL::Plane_3<R>        Plane_3;
  typedef GRAPH<SPoint_2,SSegment_2> Sphere_map;

  Sphere_map G;
  leda_node_map<leda_edge>   E;

public:

leda_sphere_map_overlayer() : G(),E(G) {}

const Sphere_map& sphere_map() const { return G; }

template <typename Iterator>
void subdivide(Iterator start, Iterator end)
/* subdivision is done in phases
   - first we partition all segments into the pieces in the
     closed postive xy-halfspace and into the pieces in the
     negative xy-halfspace
   - we sweep both halfspheres separate. Note that the boundary
     carries the same topology
   - we unify the graphs embedded into both halfspheres at
     the boundary.
*/
{
typedef leda_graph_decorator<R,Iterator> leda_graph_output;
typedef CGAL::Positive_halfsphere_geometry<R> PH_geometry;
typedef CGAL::Segment_overlay_traits<
          Iterator, leda_graph_output, PH_geometry>  PHS_traits;
typedef CGAL::generic_sweep<PHS_traits> Positive_halfsphere_sweep;

typedef CGAL::Negative_halfsphere_geometry<R> NH_geometry;
typedef CGAL::Segment_overlay_traits<
          Iterator, leda_graph_output, NH_geometry> NHS_traits;
typedef CGAL::generic_sweep<NHS_traits> Negative_halfsphere_sweep;

  std::list<SSegment_2> Lp,Lm;
  partition_xy( start, end, Lp , +1);
  partition_xy( start, end, Lm , -1);
  // both lists initialized with four quarter segments
  // supporting the xy-equator thereby separating the
  // two halfspheres
  // all other segments in the range are split into their
  // connected components with respect to the xy-plane.

  leda_node v1,v2;
  leda_graph_output O(G,E);
  typedef typename PHS_traits::INPUT Input_range;
  Positive_halfsphere_sweep SP(Input_range(Lp.begin(),Lp.end()),O);
  SP.sweep();
  CGAL_NEF_TRACEN("POS SWEEP\n"<<(dump(std::cerr),""));
  v1=G.first_node(); v2=G.last_node();
  Negative_halfsphere_sweep SM(Input_range(Lm.begin(),Lm.end()),O);
  SM.sweep();
  CGAL_NEF_TRACEN("NEG SWEEP\n"<<(dump(std::cerr),""));
  v2 = G.succ_node(v2);
  // now two CCs of sphere graph calculated
  // v1 = first node of CC in positive xy-sphere
  // v2 = first node of CC in negative xy-sphere

  merging_halfspheres(v1,v2);
  clean_trivial_sface_cycles();
  if (!Is_Plane_Map(G)) error_handler(1,"Sphere map: embedding wrong.");
  compute_faces();
}

void merge_nodes(leda_edge e1, leda_edge e2)
// e1 and e2 are two edges of the xy equator such that
// e1 is part of the positive xy-sphere bounding outer face
// e2 is part of the negative xy-sphere bounding outer face
// e1 and e2 are oppositely oriented
// the outer faces are left of the edges
// the edges are embedded orderpreserving ccw
// then the following code merges the edges of A(target(e2))
// to A(source(e1)) preserving the embedding
// afterwards source(e1) carries all edges, target(e2) is isolated
{
  leda_node v = source(e1);
  leda_edge e_pos = e1;
  leda_edge e = G.reversal(e2), er = e2;
  leda_edge e_end = e;
  do {
    leda_edge e_next = G.cyclic_adj_succ(e);
    G.move_edge(e,e_pos,target(e));
    G.move_edge(er,er,v);
    e_pos = e;
    e = e_next; er = G.reversal(e);
  } while ( e != e_end );
}

void merging_halfspheres(leda_node v1, leda_node v2)
// v1 and v2 are the definite nodes of both CCs
// where the negative y-axis pierces the sphere.
// the faces are left of edges
// edges are embedded orderpreserving ccw
{
  CGAL_NEF_TRACEN("Merging Halfspheres");
  leda_edge e1,e2,e3,e4,e1n,e2n;
  forall_sadj_edges(e1,v1)
    if ( G[target(e1)].hz()==0 && G[target(e1)].hx()<0 ) break;
  forall_sadj_edges(e2,v2)
    if ( G[target(e2)].hz()==0 && G[target(e2)].hx()>0 ) break;
  e3 = G.face_cycle_pred(e1);
  e4 = e2; e2 = G.face_cycle_pred(e2);
  while ( e1 != e3 || e2 != e4 ) {
      CGAL_NEF_TRACEN(G[source(e1)]<<" "<< G[target(e2)]);
    e1n = G.face_cycle_succ(e1); e2n = G.face_cycle_pred(e2);
    merge_nodes(e1,e2);
    e1 = e1n; e2 = e2n;
  }
}


void clean_trivial_sface_cycles()
// removes trivial face cycles at equator
// removes isolated vertices stemming from
// equator unification
{
  leda_edge_map<bool> known(G,false);
  leda_list<leda_edge> L;
  leda_list<edge_pair> Lr;
  leda_edge e;
  forall_sedges(e,G) {
    if (known[e]) continue;
    leda_edge en = G.face_cycle_succ(e);
    if ( G.face_cycle_succ(en) != e )
      continue;
    // e in trivial face cycle
    L.append(e); L.append(en);
    CGAL_NEF_TRACEN("tivial cycle "<<G[source(e)]<<G[target(e)]);
    known[e] = known[en] = true;
    leda_edge er = G.reversal(e);
    leda_edge enr = G.reversal(en);
    Lr.append(edge_pair(er,enr));
  }
  edge_pair ep;
  forall(ep,Lr) G.set_reversal(ep.first,ep.second);
  G.del_edges(L);
  leda_node v;
  forall_snodes(v,G) if ( G.outdeg(v)==0 ) G.del_node(v);
}

void compute_faces()
{
  G.compute_faces();
  leda_face f;
  leda_edge e;
  forall_sfaces(f,G) {
    CGAL_NEF_TRACEN("FACE:");
    forall_sface_edges(e,f)
      CGAL_NEF_TRACEN("  "<<SSegment_2(G[source(e)],G[target(e)]));
  }
}

void dump(std::ostream& os, leda_node v, bool nl=true) const
{  os << " ["<< ::index(v)<<"] "<<G[v];
   if (nl) os << std::endl; }

void dump(std::ostream& os) const
{
  leda_node v;
  leda_edge e;
  forall_snodes(v,G) {
    dump(os,v);
    forall_sadj_edges(e,v) {
      os << "   ->";
      dump(os,target(e),false);
      os <<" ["<<G[e]<<" ]\n";
    }
  }

}


}; // leda_sphere_map_overlayer<R>


#endif //CGAL_LEDA_SPHERE_MAP_H
