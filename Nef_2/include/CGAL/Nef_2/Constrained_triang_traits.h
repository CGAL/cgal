// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_PM_CONSTR_TRIANG_TRAITS_H
#define CGAL_PM_CONSTR_TRIANG_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/generic_sweep.h>
#include <CGAL/Nef_2/PM_checker.h>
#include <cstdlib>
#include <string>
#include <map>
#include <set>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 19
#include <CGAL/Nef_2/debug.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4355) // complaint about using 'this' to
#endif                          // initialize a member

namespace CGAL {

struct Do_nothing {
Do_nothing() {}
template <typename ARG>
void operator()(ARG&) const {}
};


template <typename PMDEC, typename GEOM, 
          typename NEWEDGE = Do_nothing>
class Constrained_triang_traits : public PMDEC {
public:
  typedef Constrained_triang_traits<PMDEC,GEOM,NEWEDGE> Self;
  typedef PMDEC                                         Base;

  // the types interfacing the sweep:
  typedef NEWEDGE                   INPUT;
  typedef typename PMDEC::Plane_map OUTPUT; 
  typedef GEOM                      GEOMETRY;

  typedef typename GEOM::Point_2     Point;
  typedef typename GEOM::Segment_2   Segment;
  typedef typename GEOM::Direction_2 Direction;

  typedef typename Base::Halfedge_handle   Halfedge_handle;
  typedef typename Base::Vertex_handle     Vertex_handle;
  typedef typename Base::Face_handle       Face_handle;
  typedef typename Base::Halfedge_iterator Halfedge_iterator;
  typedef typename Base::Vertex_iterator   Vertex_iterator;
  typedef typename Base::Face_iterator     Face_iterator;
  typedef typename Base::Halfedge_base     Halfedge_base;
  typedef typename Base::Halfedge_around_vertex_circulator
          Halfedge_around_vertex_circulator;


  using Base::point;
  using Base::is_isolated;
  using Base::first_out_edge;
  using Base::last_out_edge;
  using Base::source;
  using Base::target;
  using Base::twin;
  using Base::next;
  using Base::previous;
  using Base::cyclic_adj_succ;
  using Base::cyclic_adj_pred;
  using Base::delete_vertex;
  using Base::make_first_out_edge;

  class lt_edges_in_sweepline : public PMDEC
  {  const Point& p;
     const Halfedge_handle& e_bottom;
     const Halfedge_handle& e_top;
     const GEOMETRY& K;
    using PMDEC::point;
    using PMDEC::source;
    using PMDEC::target;
  public:
  lt_edges_in_sweepline(const Point& pi, 
     const Halfedge_handle& e1, const Halfedge_handle& e2, 
     const PMDEC& D, const GEOMETRY& k) : 
       PMDEC(D), p(pi), e_bottom(e1), e_top(e2), K(k) {}

  lt_edges_in_sweepline(const lt_edges_in_sweepline& lt) : 
     PMDEC(lt), p(lt.p), e_bottom(lt.e_bottom), e_top(lt.e_top), K(lt.K) {}

  Segment seg(const Halfedge_handle& e) const
  { return K.construct_segment(point(source(e)),point(target(e))); }

  int orientation(Halfedge_handle e, const Point& p) const
  { return K.orientation(point(source(e)),point(target(e)),p); }

  bool operator()(const Halfedge_handle& e1, const Halfedge_handle& e2) const
  { // Precondition:
    // [[p]] is identical to the source of either [[e1]] or [[e2]]. 
    if (e1 == e_bottom || e2 == e_top) return true;
    if (e2 == e_bottom || e1 == e_top) return false;
    if ( e1 == e2 ) return 0;
    int s = 0;
    if ( p == point(source(e1)) )      s =   orientation(e2,p);
    else if ( p == point(source(e2)) ) s = - orientation(e1,p);
    else CGAL_error_msg("compare error in sweep.");
    if ( s || source(e1) == target(e1) || source(e2) == target(e2) ) 
      return ( s < 0 );
    s = orientation(e2,point(target(e1)));
    if (s==0) CGAL_error_msg("parallel edges not allowed.");
    return ( s < 0 );
  }


  }; // lt_edges_in_sweepline

  class lt_pnts_xy : public PMDEC
  { const GEOMETRY& K;
  public:
    using PMDEC::point;

   lt_pnts_xy(const PMDEC& D, const GEOMETRY& k) : PMDEC(D), K(k) {}
   lt_pnts_xy(const lt_pnts_xy& lt) : PMDEC(lt), K(lt.K) {}
   bool operator()(const Vertex_handle& v1, const Vertex_handle& v2) const
   { return K.compare_xy(point(v1),point(v2)) < 0; }
  }; // lt_pnts_xy


    typedef std::map<Halfedge_handle, Halfedge_handle, lt_edges_in_sweepline> 
            Sweep_status_structure; 
    typedef typename Sweep_status_structure::iterator   ss_iterator;
    typedef typename Sweep_status_structure::value_type ss_pair;
    typedef std::set<Vertex_iterator,lt_pnts_xy> Event_Q;
    typedef typename Event_Q::const_iterator event_iterator;

    const GEOMETRY&         K;
    Event_Q                 event_Q;
    event_iterator          event_it;         
    Vertex_handle           event;
    Point                   p_sweep;
    Sweep_status_structure  SL;
    CGAL::Unique_hash_map<Halfedge_handle,ss_iterator> SLItem;
    const NEWEDGE&          Treat_new_edge;
    Halfedge_handle         e_low,e_high; // framing edges !
    Halfedge_handle         e_search;

    Constrained_triang_traits(const INPUT& in, OUTPUT& out, const GEOMETRY& k) 
      : Base(out), K(k), event_Q(lt_pnts_xy(*this,K)), 
        SL(lt_edges_in_sweepline(p_sweep,e_low,e_high,*this,K)), 
        SLItem(SL.end()),  Treat_new_edge(in)
    { CGAL_NEF_TRACEN("Constrained Triangulation Sweep"); }


  Halfedge_handle new_bi_edge(Vertex_handle v1, Vertex_handle v2)
  { // appended at v1 and v2 adj list
    Halfedge_handle e = Base::new_halfedge_pair(v1,v2);
    Treat_new_edge(e);
    return e;
  }

  Halfedge_handle new_bi_edge(Halfedge_handle e_bf, Halfedge_handle e_af)
  { // ccw before e_bf and after e_af 
    Halfedge_handle e = Base::new_halfedge_pair(e_bf,e_af,Halfedge_base(),
      Base::BEFORE, Base::AFTER);
    Treat_new_edge(e);
    return e;
  }

  Halfedge_handle new_bi_edge(Vertex_handle v, Halfedge_handle e_bf)
  { // appended at v's adj list and before e_bf
    Halfedge_handle e = Base::new_halfedge_pair(v,e_bf,Halfedge_base(),
      Base::BEFORE);
    Treat_new_edge(e);
    return e;
  }

  Segment seg(Halfedge_handle e) const
  { return K.construct_segment(point(source(e)),point(target(e))); }

  Direction dir(Halfedge_handle e) const
  { return K.construct_direction(point(source(e)),point(target(e))); }

  bool is_forward(Halfedge_handle e) const
  { return K.compare_xy(point(source(e)),point(target(e))) < 0; }




  bool edge_is_visible_from(Vertex_handle v, Halfedge_handle e)
  {
    Point p =  point(v);
    Point p1 = point(source(e));
    Point p2 = point(target(e));
    return ( K.orientation(p1,p2,p)>0 ); // left_turn
  }

  void triangulate_up(Halfedge_handle& e_apex)
  {
    CGAL_NEF_TRACEN("triangulate_up "<<seg(e_apex));
    Vertex_handle v_apex = source(e_apex);
    while (true) {
      Halfedge_handle e_vis = previous(twin(e_apex));
      bool in_sweep_line = (SLItem[e_vis] != SL.end()); 
      bool not_visible = !edge_is_visible_from(v_apex,e_vis);
        CGAL_NEF_TRACEN(" checking "<<in_sweep_line<<not_visible<<" "<<seg(e_vis));
      if ( in_sweep_line || not_visible) {
        CGAL_NEF_TRACEN("  STOP"); return;
      }
      Halfedge_handle e_back = new_bi_edge(e_apex,e_vis);
      if ( !is_forward(e_vis) ) make_first_out_edge(twin(e_back));
      e_apex = e_back;
      CGAL_NEF_TRACEN(" produced " << seg(e_apex));
    }
  }

  void triangulate_down(Halfedge_handle& e_apex)
  {
    CGAL_NEF_TRACEN("triangulate_down "<<seg(e_apex));
    Vertex_handle v_apex = source(e_apex);
    while (true) {
      Halfedge_handle e_vis = next(e_apex);
      bool in_sweep_line = (SLItem[e_vis] != SL.end());
      bool not_visible = !edge_is_visible_from(v_apex,e_vis);
        CGAL_NEF_TRACEN(" checking "<<in_sweep_line<<not_visible<<" "<<seg(e_vis));
      if ( in_sweep_line || not_visible) {
          CGAL_NEF_TRACEN("  STOP"); return;
      }
      Halfedge_handle e_vis_rev = twin(e_vis);
      Halfedge_handle e_forw = new_bi_edge(e_vis_rev,e_apex);
      e_apex = twin(e_forw);
      CGAL_NEF_TRACEN(" produced " << seg(e_apex));
    }
  }

  void triangulate_between(Halfedge_handle e_upper, Halfedge_handle e_lower)
  {
    // we triangulate the interior of the whole chain between
    // target(e_upper) and target(e_lower)
    CGAL_assertion(source(e_upper)==source(e_lower));
    CGAL_NEF_TRACE("triangulate_between\n   "<<seg(e_upper));
    CGAL_NEF_TRACEN("\n   "<<seg(e_lower));
    Halfedge_handle e_end = twin(e_lower);
    while (true) {
      Halfedge_handle e_vis =  next(e_upper);
      Halfedge_handle en_vis = next(e_vis);
      CGAL_NEF_TRACEN(" working on base e_vis " << seg(e_vis));
      CGAL_NEF_TRACEN(" next is " << seg(en_vis));
      if (en_vis == e_end) return;
      e_upper = twin(new_bi_edge(twin(e_vis),e_upper));
      CGAL_NEF_TRACEN(" produced " << seg(e_upper));
    } 
  }

  void process_event() 
  {
      CGAL_NEF_TRACEN("\nPROCESS_EVENT " << p_sweep);
    Halfedge_handle e, ep, eb_low, eb_high, e_end;
    if ( !is_isolated(event) ) {
      e = last_out_edge(event);
      ep = first_out_edge(event);
    }
    ss_iterator sit_pred, sit;
    /* PRECONDITION:
       only ingoing => e is lowest in ingoing bundle
       only outgoing => e is highest in outgoing bundle
       ingoing and outgoing => e is lowest in ingoing bundle */
    eb_high = e_end = ep;
    eb_low = e;
    CGAL_NEF_TRACEN("determining handle in SL");
    if ( e != Halfedge_handle() ) {
      point(target(e_search)) = p_sweep; // degenerate loop edge
      sit_pred = SLItem[e];
      if ( sit_pred != SL.end())  sit = --sit_pred;
      else  sit = sit_pred = --SL.upper_bound(e_search);
    } else { // event is isolated vertex
      point(target(e_search)) = p_sweep; // degenerate loop edge
      sit_pred = --SL.upper_bound(e_search);
    }

    bool ending_edges(0), starting_edges(0);
    while ( e != Halfedge_handle() ) { // walk adjacency list clockwise
      if ( SLItem[e] != SL.end() ) 
      {
        CGAL_NEF_TRACEN("ending " << seg(e));
        if (ending_edges) triangulate_between(e,cyclic_adj_succ(e));
        ending_edges = true;
        SL.erase(SLItem[e]);
        link_bi_edge_to(e,SL.end());
        // not in SL anymore
      }

      else
      {
        CGAL_NEF_TRACEN("starting "<<seg(e));
        sit = SL.insert(sit,ss_pair(e,e));
        link_bi_edge_to(e,sit);
        if ( !starting_edges ) eb_high = cyclic_adj_succ(e);
        starting_edges = true;
      }

      if (e == e_end) break;
      e = cyclic_adj_pred(e);
    }
    if (!ending_edges) 
    {
      Halfedge_handle e_vis = sit_pred->second;
      Halfedge_handle e_vis_n = cyclic_adj_succ(e_vis);
      eb_low = eb_high = new_bi_edge(event,e_vis_n); 
      CGAL_NEF_TRACEN(" producing link "<<seg(eb_low)<<"\n    before "<<seg(e_vis_n));
    }

      


    triangulate_up(eb_high);
    triangulate_down(eb_low);
    sit_pred->second = eb_low;
  }

  bool event_exists() 
  { if ( event_it != event_Q.end() ) {
      // event is set at end of loop and in init
      event = *event_it;
      p_sweep = point(event);
      return true;
    }
    return false; 
  }

  void procede_to_next_event() 
  { ++event_it; }

  void link_bi_edge_to(Halfedge_handle e, ss_iterator sit) {
    SLItem[e] = SLItem[twin(e)] = sit; 
  }

  void initialize_structures()
  {
      CGAL_NEF_TRACEN("initialize_structures ");
    
    for ( event=this->vertices_begin(); event != this->vertices_end(); ++event )
      event_Q.insert(event); // sorted order of vertices

    event_it = event_Q.begin();
    if ( event_Q.empty() ) return;
    event = *event_it;
    p_sweep = point(event); 
    if ( !is_isolated(event) ) {
      Halfedge_around_vertex_circulator 
        e(first_out_edge(event)), eend(e);
      CGAL_For_all(e,eend) {
        CGAL_NEF_TRACEN("init with "<<PE(e));
        ss_iterator sit = SL.insert(ss_pair(e,e)).first;
        link_bi_edge_to(e,sit);
      }
    }


    Vertex_handle v_tmp = this->new_vertex(); point(v_tmp) = Point();
    e_high = Base::new_halfedge_pair(event,v_tmp);
    e_low  = Base::new_halfedge_pair(event,v_tmp);
    // this are two symbolic edges just accessed as sentinels
    // they carry no geometric information
    e_search = Base::new_halfedge_pair(v_tmp,v_tmp);
    // this is just a loop used for searches in SL

    ss_iterator sit_high = SL.insert(ss_pair(e_high,e_high)).first;
    ss_iterator sit_low  = SL.insert(ss_pair(e_low,e_low)).first;
    // inserting sentinels into SL
    link_bi_edge_to(e_high, sit_high);
    link_bi_edge_to(e_low , sit_low);
    // we mark them being in the sweepline, which they will never leave 


    // we move to the second vertex:
    procede_to_next_event();
    event_exists(); // sets p_sweep for check invariants
    CGAL_NEF_TRACEN("EOF initialization");
  }

  void complete_structures() 
  {
    if (e_low != Halfedge_handle()) {
      delete_vertex(target(e_search));
    } // removing sentinels and e_search
  }


  void check_ccw_local_embedding() const
  { PM_checker<PMDEC,GEOM> C(*this,K); 
    C.check_order_preserving_embedding(event);
  }

  void check_invariants()
  {
  #ifdef CGAL_CHECK_EXPENSIVE
    if ( event_it == event_Q.end() ) return;
    check_ccw_local_embedding();
  #endif
  }

  void check_final()
  {
  #ifdef CGAL_CHECK_EXPENSIVE
    PM_checker<PMDEC,GEOM> C(*this,K); C.check_is_triangulation();
  #endif
  }

}; // Constrained_triang_traits<PMDEC,GEOM,NEWEDGE>

} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_PM_CONSTR_TRIANG_TRAITS_H
