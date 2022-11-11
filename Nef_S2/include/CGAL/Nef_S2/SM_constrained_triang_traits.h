// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_SM_CONSTRAINED_TRIANG_TRAITS_H
#define CGAL_SM_CONSTRAINED_TRIANG_TRAITS_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/generic_sweep.h>
//#include <CGAL/Nef_S2/SM_checker.h>
#include <string>
#include <map>
#include <set>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 139
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

/* For a detailed documentation see the MPI research report 2001-1-003
   which documents the planar flavor of this baby; only minor deviations
   are included in this code */

template <typename Decorator_, typename Kernel_>
class SM_constrained_triang_traits : public Decorator_ {
public:
  typedef SM_constrained_triang_traits<Decorator_,Kernel_> Self;
  typedef Decorator_                                       Base;


  typedef typename Kernel_::Point_2     Point;
  typedef typename Kernel_::Segment_2   Segment;

  typedef typename Base::SHalfedge_handle   SHalfedge_handle;
  typedef typename Base::SVertex_handle     SVertex_handle;
  typedef typename Base::SFace_handle       SFace_handle;
  typedef typename Base::SHalfedge_iterator SHalfedge_iterator;
  typedef typename Base::SVertex_iterator   SVertex_iterator;
  typedef typename Base::SFace_iterator     SFace_iterator;
  typedef typename Base::SHalfedge_around_svertex_circulator
          SHalfedge_around_svertex_circulator;

  // the types interfacing the sweep:
  typedef std::pair<SVertex_iterator,SVertex_iterator> INPUT;
  typedef Decorator_                                 OUTPUT;
  typedef Kernel_                                    GEOMETRY;

  using Base::is_isolated;
  using Base::first_out_edge;
  using Base::last_out_edge;
  using Base::new_svertex;
  using Base::assoc_info;
  using Base::incident_mark;
  using Base::cyclic_adj_succ;
  using Base::cyclic_adj_pred;
  using Base::delete_vertex;

  class lt_edges_in_sweepline : public Decorator_
  {  const Point& p;
     const SHalfedge_handle& e_bottom;
     const SHalfedge_handle& e_top;
     const Kernel_& K;
  public:
  lt_edges_in_sweepline(const Point& pi,
     const SHalfedge_handle& e1, const SHalfedge_handle& e2,
     const Decorator_& D, const Kernel_& k) :
       Decorator_(D), p(pi), e_bottom(e1), e_top(e2), K(k) {}

  lt_edges_in_sweepline(const lt_edges_in_sweepline& lt) :
     Decorator_(lt), p(lt.p),
     e_bottom(lt.e_bottom), e_top(lt.e_top), K(lt.K) {}

  Segment seg(const SHalfedge_handle& e) const
  { return K.construct_segment(point(source(e)),point(target(e))); }

  int orientation(SHalfedge_handle e, const Point& p) const
  { return K.orientation(e->source()->point(),e->twin()->source()->point(),p); }

  bool operator()(const SHalfedge_handle& e1, const SHalfedge_handle& e2) const
  { // Precondition:
    // [[p]] is identical to the source of either [[e1]] or [[e2]].
    if (e1 == e_bottom || e2 == e_top) return true;
    if (e2 == e_bottom || e1 == e_top) return false;
    if ( e1 == e2 ) return 0;
    int s = 0;
    if ( p == e1->source()->point() )
      s =   orientation(e2,p);
    else if ( p == e2->source()->point() )
      s = - orientation(e1,p);
    else CGAL_error_msg("compare error in sweep.");
    if ( s || e1->source() == e1->twin()->source() ||
         e2->source() == e2->twin()->source())
      return ( s < 0 );
    s = orientation(e2,e1->twin()->source()->point());
    if (s==0) CGAL_error_msg("parallel edges not allowed.");
    return ( s < 0 );
  }


  }; // lt_edges_in_sweepline

  class lt_pnts_xy : public Decorator_
  { const Kernel_& K;
  public:
   lt_pnts_xy(const Decorator_& D, const Kernel_& k) : Decorator_(D), K(k) {}
   lt_pnts_xy(const lt_pnts_xy& lt) : Decorator_(lt), K(lt.K) {}
   int operator()(const SVertex_handle& v1, const SVertex_handle& v2) const
   { return K.compare_xy(v1->point(),v2->point()) < 0; }
  }; // lt_pnts_xy


  typedef std::map<SHalfedge_handle, SHalfedge_handle, lt_edges_in_sweepline>
          Sweep_status_structure;
  typedef typename Sweep_status_structure::iterator   ss_iterator;
  typedef typename Sweep_status_structure::value_type ss_pair;
  typedef std::set<SVertex_iterator,lt_pnts_xy> Event_Q;
  typedef typename Event_Q::const_iterator event_iterator;

  const GEOMETRY&         K;
  Event_Q                 event_Q;
  event_iterator          event_it;
  SVertex_handle           event;
  Point                   p_sweep;
  Sweep_status_structure  SL;
  CGAL::Unique_hash_map<SHalfedge_handle,ss_iterator> SLItem;
  SHalfedge_handle         e_low,e_high; // framing edges !
  SHalfedge_handle         e_search;
  SVertex_iterator         v_first, v_beyond;

  SM_constrained_triang_traits(const INPUT& in, OUTPUT& out,
                               const GEOMETRY& k)
      : Base(out), K(k), event_Q(lt_pnts_xy(*this,K)),
        SL(lt_edges_in_sweepline(p_sweep,e_low,e_high,*this,K)),
        SLItem(SL.end()), v_first(in.first), v_beyond(in.second)
    { CGAL_NEF_TRACEN("Constrained Triangulation Sweep"); }

  /* |treat_new_sedge| is used to forward information that exists
     at input edges of the triangulation as such it spreads input
     information to the newly created edges of the triangulation;
     the used operation incident_mark refers to the base class of
     |*this| */
  void treat_new_edge(SHalfedge_handle e)
  { assoc_info(e);
    e->mark() = incident_mark(e) = incident_mark(e->twin()) =
      incident_mark(e->snext());
    CGAL_NEF_TRACEN(" treat_new_edge "<<PH(e));
  }

  SHalfedge_handle new_bi_edge(SVertex_handle v1, SVertex_handle v2)
  { // appended at v1 and v2 adj list
    SHalfedge_handle e = Base::new_shalfedge_pair(v1,v2);
    treat_new_edge(e);
    return e;
  }

  SHalfedge_handle new_bi_edge(SHalfedge_handle e_bf, SHalfedge_handle e_af)
  { // ccw before e_bf and after e_af
    SHalfedge_handle e =
      Base::new_shalfedge_pair(e_bf,e_af,Base::BEFORE, Base::AFTER);
    treat_new_edge(e);
    return e;
  }

  SHalfedge_handle new_bi_edge(SVertex_handle v, SHalfedge_handle e_bf)
  { // appended at v's adj list and before e_bf
    SHalfedge_handle e = Base::new_shalfedge_pair(v,e_bf,Base::BEFORE);
    treat_new_edge(e);
    return e;
  }

  Segment segment(SHalfedge_handle e) const
  { return K.construct_segment(e->source()->point(),e->twin()->source()->point()); }

  bool is_forward(SHalfedge_handle e) const
  { return K.compare_xy(point(source(e)),point(target(e))) < 0; }

  bool edge_is_visible_from(SVertex_handle v, SHalfedge_handle e)
  {
    Point p =  v->point();
    Point p1 = e->source()->point();
    Point p2 = e->twin()->source()->point();
    return ( K.orientation(p1,p2,p) > 0 ); // left_turn
  }

  void triangulate_up(SHalfedge_handle& e_apex)
  {
    CGAL_NEF_TRACEN("triangulate_up "<<segment(e_apex));
    SVertex_handle v_apex = e_apex->source();
    while (true) {
      SHalfedge_handle e_vis = e_apex->twin()->sprev();
      bool in_sweep_line = (SLItem[e_vis] != SL.end());
      bool not_visible = !edge_is_visible_from(v_apex,e_vis);
        CGAL_NEF_TRACEN(" checking "<<in_sweep_line<<not_visible<<" "<<segment(e_vis));
      if ( in_sweep_line || not_visible) {
        CGAL_NEF_TRACEN("  STOP"); return;
      }
      SHalfedge_handle e_back = new_bi_edge(e_apex,e_vis);
      e_apex = e_back;
      CGAL_NEF_TRACEN(" produced " << segment(e_apex));
    }
  }

  void triangulate_down(SHalfedge_handle& e_apex)
  {
    CGAL_NEF_TRACEN("triangulate_down "<<segment(e_apex));
    SVertex_handle v_apex = e_apex->source();
    while (true) {
      SHalfedge_handle e_vis = e_apex->snext();
      bool in_sweep_line = (SLItem[e_vis] != SL.end());
      bool not_visible = !edge_is_visible_from(v_apex,e_vis);
        CGAL_NEF_TRACEN(" checking "<<in_sweep_line<<not_visible<<" "<<segment(e_vis));
      if ( in_sweep_line || not_visible) {
          CGAL_NEF_TRACEN("  STOP"); return;
      }
      SHalfedge_handle e_vis_rev = e_vis->twin();
      SHalfedge_handle e_forw = new_bi_edge(e_vis_rev,e_apex);
      e_apex = e_forw->twin();
      CGAL_NEF_TRACEN(" produced " << segment(e_apex));
    }
  }

  void triangulate_between(SHalfedge_handle e_upper, SHalfedge_handle e_lower)
  {
    // we triangulate the interior of the whole chain between
    // target(e_upper) and target(e_lower)
    CGAL_assertion(e_upper->source()==e_lower->source());
    CGAL_NEF_TRACE("triangulate_between\n   "<<segment(e_upper));
    CGAL_NEF_TRACEN("\n   "<<segment(e_lower));
    SHalfedge_handle e_end = e_lower->twin();
    while (true) {
      SHalfedge_handle e_vis =  e_upper->snext();
      SHalfedge_handle en_vis = e_vis->snext();
      CGAL_NEF_TRACEN(" working on base e_vis " << segment(e_vis));
      CGAL_NEF_TRACEN(" next is " << segment(en_vis));
      if (en_vis == e_end) return;
      e_upper = new_bi_edge(e_vis->twin(),e_upper)->twin();
      CGAL_NEF_TRACEN(" produced " << segment(e_upper));
    }
  }

  void process_event()
  {
      CGAL_NEF_TRACEN("\nPROCESS_EVENT " << p_sweep);
    SHalfedge_handle e, ep, eb_low, eb_high, e_end;
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
    if ( e != SHalfedge_handle() ) {
      e_search->twin()->source()->point() = p_sweep; // degenerate loop edge
      sit_pred = SLItem[e];
      if ( sit_pred != SL.end())  sit = --sit_pred;
      else  sit = sit_pred = --SL.upper_bound(e_search);
    } else { // event is isolated vertex
      e_search->twin()->source()->point() = p_sweep; // degenerate loop edge
      sit_pred = --SL.upper_bound(e_search);
    }

    bool ending_edges(0), starting_edges(0);
    while ( e != SHalfedge_handle() ) { // walk adjacency list clockwise
      if ( SLItem[e] != SL.end() )
      {
        CGAL_NEF_TRACEN("ending " << segment(e));
        if (ending_edges) triangulate_between(e,cyclic_adj_succ(e));
        ending_edges = true;
        SL.erase(SLItem[e]);
        link_bi_edge_to(e,SL.end());
        // not in SL anymore
      }
      else
      {
        CGAL_NEF_TRACEN("starting "<<segment(e));
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
      SHalfedge_handle e_vis = sit_pred->second;
      SHalfedge_handle e_vis_n = cyclic_adj_succ(e_vis);
      eb_low = eb_high = new_bi_edge(event,e_vis_n);
      CGAL_NEF_TRACEN(" producing link "<<segment(eb_low)<<
             "\n    before "<<segment(e_vis_n));
    }

    triangulate_up(eb_high);
    triangulate_down(eb_low);
    sit_pred->second = eb_low;
  }

  bool event_exists()
  { if ( event_it != event_Q.end() ) {
      // event is set at end of loop and in init
      event = *event_it;
      p_sweep = event->point();
      return true;
    }
    return false;
  }

  void procede_to_next_event()
  { ++event_it; }

  void link_bi_edge_to(SHalfedge_handle e, ss_iterator sit) {
    SLItem[e] = SLItem[e->twin()] = sit;
  }

  void initialize_structures()
  {
      CGAL_NEF_TRACEN("initialize_structures ");

    for ( event = v_first; event != v_beyond; ++event )
      event_Q.insert(event); // sorted order of vertices

    event_it = event_Q.begin();
    if ( event_Q.empty() ) return;
    event = *event_it;
    p_sweep = event->point();
    if ( !is_isolated(event) ) {
      SHalfedge_around_svertex_circulator
        e(first_out_edge(event)), eend(e);
      CGAL_For_all(e,eend) {
        CGAL_NEF_TRACEN("init with "<<PH(e));
        ss_iterator sit = SL.insert(ss_pair(e,e)).first;
        link_bi_edge_to(e,sit);
      }
    }


    SVertex_handle v_tmp = new_svertex(Point());
    e_high = Base::new_shalfedge_pair(event,v_tmp);
    e_low  = Base::new_shalfedge_pair(event,v_tmp);
    // this are two symbolic edges just accessed as sentinels
    // they carry no geometric information
    e_search = Base::new_shalfedge_pair(v_tmp,v_tmp);
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
    if (e_low != SHalfedge_handle()) {
      delete_vertex(e_search->twin()->source());
    } // removing sentinels and e_search
  }


  void check_ccw_local_embedding() const
  {
#ifdef CGAL_CHECK_EXPENSIVEXXX
    PM_checker<Decorator_,Kernel_> C(*this,K);
    C.check_order_preserving_embedding(event);
#endif
  }

  void check_invariants()
  {
#ifdef CGAL_CHECK_EXPENSIVEXXX
    if ( event_it == event_Q.end() ) return;
    check_ccw_local_embedding();
#endif
  }

  void check_final()
  {
#ifdef CGAL_CHECK_EXPENSIVEXXX
    PM_checker<Decorator_,Kernel_> C(*this,K); C.check_is_triangulation();
#endif
  }

}; // SM_constrained_triang_traits<Decorator_,Kernel_,New_edge_>


} //namespace CGAL
#endif // CGAL_SM_CONSTRAINED_TRIANG_TRAITS_H
