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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_SEGMENT_OVERLAY_TRAITS_H
#define CGAL_SEGMENT_OVERLAY_TRAITS_H

#include <CGAL/license/Nef_2.h>


#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 23
#include <CGAL/Nef_2/debug.h>

#if defined(CGAL_USE_LEDA_LIBRARY)
#include <CGAL/LEDA_basic.h>

#include <LEDA/core/tuple.h>
#include <LEDA/core/slist.h>
#include <LEDA/core/list.h>
#include <LEDA/core/map.h>
#include <LEDA/core/map2.h>
#include <LEDA/core/sortseq.h>
#include <LEDA/core/p_queue.h>
#include <LEDA/core/impl/ab_tree.h>
#include <LEDA/core/impl/bb_tree.h>
#include <LEDA/core/impl/rb_tree.h>
#include <LEDA/core/impl/rs_tree.h>
#include <LEDA/core/impl/skiplist.h>

#include <utility>
#include <sstream>

namespace CGAL {
#ifdef CGAL_NEF_DEBUG
#define PIS(s) (s->first())
#endif

template <typename IT, typename PMDEC, typename GEOM>
class leda_seg_overlay_traits {
public:
  typedef IT                               ITERATOR;
  typedef std::pair<IT,IT>                 INPUT;
  typedef PMDEC                            OUTPUT;
  typedef typename PMDEC::Vertex_handle    Vertex_handle;
  typedef typename PMDEC::Halfedge_handle  Halfedge_handle;
  typedef GEOM                             GEOMETRY;
  typedef typename GEOMETRY::Point_2       Point_2;
  typedef typename GEOMETRY::Segment_2     Segment_2;

  typedef leda_two_tuple<Segment_2,ITERATOR> seg_pair;
  typedef seg_pair*                          ISegment;
  typedef leda_list<seg_pair>                IList;
  typedef typename IList::iterator           ilist_iterator;


  // types interfacing the generic sweep frame:
  ITERATOR its, ite;
  OUTPUT&  GO;
  const GEOMETRY& K;

  class cmp_segs_at_sweepline : public CGAL_LEDA_SCOPE::leda_cmp_base<ISegment>
  { const Point_2& p;
    ISegment s_bottom, s_top; // sentinel segments
    const GEOMETRY& K;
  public:
   cmp_segs_at_sweepline(const Point_2& pi, 
     ISegment s1, ISegment s2, const GEOMETRY& k) : 
     p(pi), s_bottom(s1), s_top(s2), K(k) {}

   int operator()(const ISegment& is1, const ISegment& is2) const
   { // Precondition: p is identical to the left endpoint of s1 or s2. 
     if ( is2 == s_top || is1 == s_bottom ) return -1;
     if ( is1 == s_top || is2 == s_bottom ) return 1;
     if ( is1 == is2 ) return 0;
     const Segment_2& s1 = is1->first();
     const Segment_2& s2 = is2->first();
     int s = 0;
     if ( p == K.source(s1) )      s =   K.orientation(s2,p);
     else if ( p == K.source(s2) ) s = - K.orientation(s1,p);
     else CGAL_error_msg("compare error in sweep.");
     if ( s || K.is_degenerate(s1) || K.is_degenerate(s2) ) 
       return s;
    
     s = K.orientation(s2,K.target(s1));
     if (s==0) return static_cast<int>( is1 - is2 );
     // overlapping segments are not equal
     return s;
   }
  };

  struct cmp_pnts_xy : public CGAL_LEDA_SCOPE::leda_cmp_base<Point_2>
  { const GEOMETRY& K;
  public:
   cmp_pnts_xy(const GEOMETRY& k) : K(k) {}
   int operator()(const Point_2& p1, const Point_2& p2) const
   { return K.compare_xy(p1,p2); }
  };

  //  typedef CGAL_LEDA_SCOPE::skiplist                              SearchTree;
  //  typedef typename SearchTree::item                              ST_item;
  typedef CGAL_LEDA_SCOPE::seq_item                         ST_item;
  typedef leda_sortseq<Point_2, ST_item>             EventQueue; 
  typedef leda_sortseq<ISegment, ST_item>            SweepStatus;
  typedef leda_p_queue<Point_2,ISegment>                         SegQueue; 
  typedef leda_map<ST_item,Halfedge_handle>    AssocEdgeMap;
  typedef leda_slist<ITERATOR>                                   IsoList;
  typedef typename IsoList::item                        slist_item;
  typedef leda_map<ST_item, IsoList* >         AssocIsoMap;
  typedef leda_map2<ISegment,ISegment,ST_item> EventHash;

    
    ST_item  event;
    Point_2                    p_sweep;
    cmp_pnts_xy                cmp;
    EventQueue                 XS;
    seg_pair                   sl,sh;
    cmp_segs_at_sweepline      SLcmp;
    SweepStatus                YS;
    SegQueue                   SQ;

    EventHash                  IEvent;
    IList                      Internal;
    AssocEdgeMap               Edge_of;
    AssocIsoMap                Isos_of;

    leda_seg_overlay_traits(const INPUT& in, OUTPUT& G, 
      const GEOMETRY& k) : 
      its(in.first), ite(in.second), GO(G), K(k),
      cmp(K), XS(cmp), SLcmp(p_sweep,&sl,&sh,K), YS(SLcmp), SQ(cmp),
      IEvent(0), Edge_of(0), Isos_of(0) {}


  leda_string dump_structures() const
  { 
    std::ostringstream out;
    out << "SQ= ";
    CGAL_LEDA_SCOPE::pq_item pqit;
    forall_items(pqit,SQ) {
      if (SQ.prio(pqit)==XS.key(XS.succ(XS.min_item()))) 
      { out << SQ.inf(pqit)->first(); }
      pqit = SQ.next_item(pqit);
    }
    ST_item sit;
    out << "\nXS=\n";
    forall_items(sit,XS)
      out << "  " << XS.key(sit) << " " << XS.inf(sit) 
          <<std::endl;
    out << "YS=\n";
    for( sit = YS.max_item(); sit; sit=YS.pred(sit) )
      out << "  "<<YS.key(sit)->first()<<" "<<YS.inf(sit)<<std::endl;
    leda_string res(out.str().c_str()); 
    return res;
  }


  Point_2 source(ISegment is) const
  { return K.source(is->first()); }
  Point_2 target(ISegment is) const
  { return K.target(is->first()); }
  ITERATOR original(ISegment s) const
  { return s->second(); }

  int orientation(ST_item sit, const Point_2& p) const
  { return K.orientation(YS.key(sit)->first(),p); }

  bool collinear(ST_item sit1, 
                 ST_item sit2) const
  { Point_2 ps = source(YS.key(sit2)), pt = target(YS.key(sit2));
    return ( orientation(sit1,ps)==0 &&
             orientation(sit1,pt)==0 );
  }


  void compute_intersection(ST_item sit0)
  {    
    ST_item sit1 = YS.succ(sit0);
    if ( sit0 == YS.min_item() || sit1 == YS.max_item() ) return;
    ISegment s0 = YS.key(sit0);
    ISegment s1 = YS.key(sit1);
    int or0 = K.orientation(s0->first(),target(s1));
    int or1 = K.orientation(s1->first(),target(s0));
    if ( or0 <= 0 && or1 >= 0  ) { 
      ST_item it = IEvent(YS.key(sit0),YS.key(sit1));
      if ( it==0 ) {
        Point_2 q = K.intersection(s0->first(),s1->first());
        it = XS.insert(q,sit0);
      }
      YS.change_inf(sit0, it);
    }
  } 

  void initialize_structures()
  {
    CGAL_NEF_TRACEN("initialize_structures");
    ITERATOR it_s;  
    for ( it_s=its; it_s != ite; ++it_s ) {
      Segment_2 s = *it_s;
      ST_item it1 = 
            XS.insert( K.source(s), ST_item(nil));
      ST_item it2 = 
            XS.insert( K.target(s), ST_item(nil));
      if (it1 == it2) {
        if ( Isos_of[it1] == 0 ) Isos_of[it1] = new IsoList;
        Isos_of[it1]->push(it_s);
        continue;  // ignore zero-length segments in SQ/YS
      }
          
      Point_2 p = XS.key(it1);
      Point_2 q = XS.key(it2);
          
      Segment_2 s1; 
      if ( K.compare_xy(p,q) < 0 ) 
        s1 = K.construct_segment(p,q);
      else
        s1 = K.construct_segment(q,p);
          
      Internal.append(seg_pair(s1,it_s));
      SQ.insert(K.source(s1),&Internal[Internal.last()]);
    }

    // insert a lower and an upper sentinel segment
    YS.insert(&sl,ST_item(nil));
    YS.insert(&sh,ST_item(nil));
    CGAL_NEF_TRACEN("end of initialization\n"<<YS.size());
  }

  bool event_exists() 
  { 
    if (!XS.empty()) { 
      // event is set at end of loop and in init
      event = (XS.min)();
      p_sweep = XS.key(event);
      return true;
    }
    return false; 
  }

  void procede_to_next_event() 
  { XS.del_item(event); }

  void process_event() 
  {
    CGAL_NEF_TRACEN("\n\n >>> process_event: "<<p_sweep<<" "<<XS[event]<<" "<<event);

    Vertex_handle v = GO.new_vertex(p_sweep);
    ST_item sit = XS.inf(event);
        
      ST_item sit_succ(0), sit_pred(0), sit_pred_succ(0), sit_first(0);
      if (sit == nil) 
        {
          Segment_2 s_sweep = K.construct_segment(p_sweep,p_sweep);
          seg_pair sp(s_sweep, ite);
          sit_succ = YS.locate( &sp );
          if ( sit_succ != YS.max_item() && 
               orientation(sit_succ,p_sweep) == 0 ) 
            sit = sit_succ;
          else  {
            sit_pred = YS.pred(sit_succ);
            sit_pred_succ = sit_succ;
          }
          CGAL_NEF_TRACEN("looked up p_sweep "<<PIS(YS.key(sit_succ)));
        }



      /* If no segment contains p_sweep then sit_pred and sit_succ are
         correctly set after the above locate operation, if a segment
         contains p_sweep sit_pred and sit_succ are set below when
         determining the bundle.*/

      if (sit != nil) { // key(sit) is an ending or passing segment
        CGAL_NEF_TRACEN("ending/passing segs");
        while ( YS.inf(sit) == event ||
                YS.inf(sit) == YS.succ(sit) ) // overlapping
          sit = YS.succ(sit);
        sit_succ = YS.succ(sit); 
        ST_item sit_last = sit;

        ST_item xit = YS.inf(sit_last);
        if (xit) { 
          ISegment s1 = YS.key(sit_last);
          ISegment s2 = YS.key(sit_succ);
          IEvent(s1,s2) = xit;
            CGAL_NEF_TRACEN("hashing "<<PIS(s1)<<PIS(s2)<<xit);
        } 
          
        bool overlapping;
        do {
          ISegment s = YS.key(sit);
          ST_item sit_next = YS.pred(sit);
          overlapping = (YS.inf(sit_next) == sit);
          Halfedge_handle e = Edge_of[sit];
          if ( !overlapping ) {
            CGAL_NEF_TRACEN("connecting edge to node "<<PIS(s)<<" "<<sit);
            GO.link_as_target_and_append(v,e);
          }
          GO.supporting_segment(e,original(s));
          if ( target(s) == p_sweep ) { // ending segment
              CGAL_NEF_TRACEN("ending segment "<<PIS(s));
            if ( overlapping ) 
	      YS.change_inf(sit_next,YS.inf(sit));
            YS.del_item(sit);
            GO.ending_segment(v,original(s));
          } else {  // passing segment
	    CGAL_NEF_TRACEN("passing segment "<<PIS(s));
            if ( YS.inf(sit) != YS.succ(sit) ) 
              YS.change_inf(sit, ST_item(0));
            GO.passing_segment(v,original(s));
          }
          sit = sit_next;
        } 
        while ( YS.inf(sit) == event || overlapping ||
                YS.inf(sit) == YS.succ(sit) );
                  
        sit_pred = sit;
        sit_first = sit_pred_succ = YS.succ(sit_pred); // first item of bundle

        CGAL_NEF_TRACE("event bundles between\n   "<<PIS(YS.key(sit_succ)));
        CGAL_NEF_TRACEN("\n   "<<PIS(YS.key(sit_pred)));

        while ( sit != sit_succ ) {
          ST_item sub_first = sit;
          ST_item sub_last  = sub_first;
                            
          while (YS.inf(sub_last) == YS.succ(sub_last))
            sub_last = YS.succ(sub_last);
                            
          if (sub_last != sub_first)
            YS.reverse_items(sub_first, sub_last);
                            
          sit = YS.succ(sub_first);
        }
                        
        // reverse the entire bundle
        if (sit_first != sit_succ) 
          YS.reverse_items(YS.succ(sit_pred),YS.pred(sit_succ));


      } // if (sit != nil)

    CGAL_assertion(sit_pred);
    GO.halfedge_below(v,Edge_of[sit_pred]);
    if ( Isos_of[event] != 0 ) {
      const IsoList& IL = *(Isos_of[event]);
      slist_item iso_it;
      for (iso_it = IL.first(); iso_it; iso_it=IL.succ(iso_it) ) 
        GO.trivial_segment(v,IL[iso_it] );
      delete (Isos_of[event]); // clean up the list
    }


    ISegment next_seg;
    CGAL_LEDA_SCOPE::pq_item next_it = SQ.find_min();
    while ( next_it && 
            (next_seg = SQ.inf(next_it), p_sweep == source(next_seg)) ) {
      ST_item s_sit = YS.locate_succ(next_seg);
      ST_item p_sit = YS.pred(s_sit);

      CGAL_NEF_TRACEN("inserting "<<PIS(next_seg)<<" at "<<PIS(YS.key(s_sit))); 
      if ( YS.max_item() != s_sit &&
           orientation(s_sit, source(next_seg) ) == 0 &&
           orientation(s_sit, target(next_seg) ) == 0 )
        sit = YS.insert_at(s_sit, next_seg, s_sit);
      else 
        sit = YS.insert_at(s_sit, next_seg, ST_item(nil));
      CGAL_assertion(YS.succ(sit)==s_sit);

      if ( YS.min_item() != p_sit &&
           orientation(p_sit, source(next_seg) ) == 0 &&
           orientation(p_sit, target(next_seg) ) == 0 )
        YS.change_inf(p_sit, sit);
      CGAL_assertion(YS.succ(p_sit)==sit);
                 
      XS.insert(target(next_seg), sit);
      GO.starting_segment(v,original(next_seg));
                 
      // delete minimum and assign new minimum to next_seg
      SQ.del_min();    
      next_it = SQ.find_min();
    }

    for( ST_item sitl = YS.pred(sit_succ); sitl != sit_pred; 
         sitl = YS.pred(sitl) ) {
      if ( YS.inf(sitl) != YS.succ(sitl) ) { // non-overlapping
        CGAL_NEF_TRACEN("non-overlapping "<<PIS(YS.key(sitl))<<" "<<sitl);
        Edge_of[sitl] = GO.new_halfedge_pair_at_source(v);
      } else {
        CGAL_NEF_TRACEN("overlapping "<<PIS(YS.key(sitl)));
        Edge_of[sitl] = Edge_of[ YS.succ(sitl) ];
      }
    }
    sit_first = YS.succ(sit_pred);


    CGAL_assertion(sit_pred); CGAL_assertion(sit_pred_succ);
    ST_item xit = YS.inf(sit_pred);
    if ( xit ) { 
      ISegment s1 = YS.key(sit_pred);
      ISegment s2 = YS.key(sit_pred_succ);
      IEvent(s1,s2) = xit;
        CGAL_NEF_TRACEN("hashing "<<PIS(s1)<<PIS(s2)<<xit);
      YS.change_inf(sit_pred, ST_item(0));
    }
          
    compute_intersection(sit_pred); 
    sit = YS.pred(sit_succ);
    if (sit != sit_pred)
      compute_intersection(sit);


  }

  void complete_structures() {}
  void check_invariants() {CGAL_NEF_TRACEN("check_invariants\n"<<dump_structures());}
  void check_final() {}

}; // leda_seg_overlay_traits

} // namespace CGAL

#endif // defined(CGAL_USE_LEDA_LIBRARY)
#if !defined(CGAL_USE_LEDA_LIBRARY)
#include <list>
#include <string>
#include <sstream>
#include <map>
#include <CGAL/Multiset.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

template <typename IT, typename PMDEC, typename GEOM>
class stl_seg_overlay_traits {
public:
  typedef IT                               ITERATOR;
  typedef std::pair<IT,IT>                 INPUT;
  typedef PMDEC                            OUTPUT;
  typedef typename PMDEC::Vertex_handle    Vertex_handle;
  typedef typename PMDEC::Halfedge_handle  Halfedge_handle;
  typedef GEOM                             GEOMETRY;
  typedef typename GEOMETRY::Point_2       Point_2;
  typedef typename GEOMETRY::Segment_2     Segment_2;

  typedef std::pair<Segment_2,ITERATOR>     seg_pair;
  typedef seg_pair*                         ISegment;
  typedef std::list<seg_pair>               IList;
  typedef typename IList::const_iterator    ilist_iterator;


  // types interfacing the generic sweep frame
  ITERATOR its, ite;
  OUTPUT&  GO;
  const GEOMETRY& K;


  class lt_segs_at_sweepline 
  { 
    const Point_2& p;
    ISegment s_bottom, s_top; // sentinel segments
    const GEOMETRY& K;

  public:
    lt_segs_at_sweepline(const Point_2& pi, 
			 ISegment s1, ISegment s2, 
			 const GEOMETRY& k) 
     : p(pi), s_bottom(s1), s_top(s2), K(k) 
    {}
    
    lt_segs_at_sweepline(const lt_segs_at_sweepline& lt) 
      : p(lt.p), s_bottom(lt.s_bottom), s_top(lt.s_top), K(lt.K) 
    {}
    
    template <typename ss_pair>
    bool
    operator()(const ISegment& is1, const ss_pair& ss2) const
    {
      return operator()(is1, ss2.first);
    }

    bool 
    operator()(const ISegment& is1, const ISegment& is2) const
    { 
      if ( is2 == s_top || is1 == s_bottom ) return true;
      if ( is1 == s_top || is2 == s_bottom ) return false;
      if ( is1 == is2 ) return false;
      // Precondition: p is contained in s1 or s2. 
      const Segment_2& s1 = is1->first;
      const Segment_2& s2 = is2->first;

      CGAL_assertion_msg(( K.orientation(s1,p) == 0 ) ||  ( K.orientation(s2,p) == 0 ) ,"compare error in sweep.");

      int s = 0;
      if( p == K.source(s1) )
	s = K.orientation(s2,p);
      else
	s = - K.orientation(s1,p);
      if ( s || K.is_degenerate(s1) || K.is_degenerate(s2) ) 
	return ( s < 0 );
      
      s = K.orientation(s2,K.target(s1));
      if (s==0) return ( is1 - is2 ) < 0;
      // overlapping segments are not equal
      return ( s < 0 );
    }
  };
  

  class compare_segs_at_sweepline 
  { 
    const Point_2& p;
    ISegment s_bottom, s_top; // sentinel segments
    const GEOMETRY& K;

    CGAL::Comparison_result o2c(int o) const {
      if(o > 0) return CGAL::LARGER;
      if(o < 0) return CGAL::SMALLER;
      return CGAL::EQUAL;
    }

  public:
    compare_segs_at_sweepline(const Point_2& pi, 
			 ISegment s1, ISegment s2, 
			      const GEOMETRY& k) 
     : p(pi), s_bottom(s1), s_top(s2), K(k) 
    {}
    
    compare_segs_at_sweepline(const compare_segs_at_sweepline& lt) 
      : p(lt.p), s_bottom(lt.s_bottom), s_top(lt.s_top), K(lt.K) 
    {}
    
    template <typename ss_pair>
    bool
    operator()(const ISegment& is1, const ss_pair& ss2) const
    {
      return operator()(is1, ss2.first);
    }

    CGAL::Comparison_result
    operator()(const ISegment& is1, const ISegment& is2) const
    { 
      if ( is2 == s_top || is1 == s_bottom ) return CGAL::SMALLER;
      if ( is1 == s_top || is2 == s_bottom ) return CGAL::LARGER;
      if ( is1 == is2 ) return CGAL::EQUAL;
      // Precondition: p is contained in s1 or s2. 
      const Segment_2& s1 = is1->first;
      const Segment_2& s2 = is2->first;

      CGAL_assertion_msg(( K.orientation(s1,p) == 0 ) ||  ( K.orientation(s2,p) == 0 ) ,"compare error in sweep.");

      int s = 0;
      if(K.is_degenerate(s1))
	return o2c(K.orientation(s2,p));
      if(K.is_degenerate(s2))
	return o2c(-K.orientation(s1,p));
      
      s = - K.orientation(s1,p);
      if(s!=0)
	return o2c(s);
      s = K.orientation(s2,p);
      if(s!=0)
	return o2c(s);
      return o2c(K.orientation(s2,K.target(s1)));
      /*
      if( p == K.source(s1) )
	s = K.orientation(s2,p);
      else
	s = - K.orientation(s1,p);
      if ( s || K.is_degenerate(s1) || K.is_degenerate(s2) )
	if(s < 0) return CGAL::SMALLER;
	else if(s > 0) return CGAL::LARGER;
	else return CGAL::EQUAL;
      
      s = K.orientation(s2,K.target(s1));
      //	if (s==0) {
      //	if(is1 < is2) return CGAL::SMALLER;
      //	if (is1 > is2) return CGAL::LARGER;
      //	return CGAL::EQUAL;
      //      }

      // overlapping segments are not equal
      if(s < 0) return CGAL::SMALLER;
      if(s > 0) return CGAL::LARGER;
      return CGAL::EQUAL;
*/
    }
  };

  class compare_pnts_xy {
    const GEOMETRY& K;
    
    public:
    compare_pnts_xy(const GEOMETRY& k) 
      : K(k) 
    {}
    
    compare_pnts_xy(const compare_pnts_xy& lt) 
      : K(lt.K) 
    {}
    
    CGAL::Comparison_result
    operator()(const Point_2& p1, const Point_2& p2) const
    { 
      int c = K.compare_xy(p1,p2);
      if(c < 0) return CGAL::SMALLER;
      if(c > 0) return CGAL::LARGER;
      return CGAL::EQUAL;
    }
  };
  
  struct lt_pnts_xy { 
    const GEOMETRY& K;
    
  public:
    lt_pnts_xy(const GEOMETRY& k) 
      : K(k) 
    {}
    
    lt_pnts_xy(const lt_pnts_xy& lt) 
      : K(lt.K) 
    {}
    
    bool 
    operator()(const Point_2& p1, const Point_2& p2) const
    { 
      return K.compare_xy(p1,p2) < 0; 
    }
  };
  
  typedef std::list<ITERATOR>                            IsoList;
  typedef CGAL::Multiset<Point_2, compare_pnts_xy>       EventQueue;
  typedef typename EventQueue::iterator                  event_iterator;
  typedef Unique_hash_map<Point_2*, IsoList*>            X2iso;

  typedef CGAL::Multiset<ISegment, compare_segs_at_sweepline>  SweepStatus;
  typedef typename SweepStatus::iterator                       ss_iterator;
  typedef typename SweepStatus::const_iterator                 ss_const_iterator;
  typedef CGAL::Unique_hash_map<ISegment, event_iterator>      Y2X;
  typedef CGAL::Unique_hash_map<ISegment, ISegment>            Y2Y;
  typedef CGAL::Unique_hash_map<Point_2*, ss_iterator>         X2Y;
  typedef CGAL::Unique_hash_map<ISegment, Halfedge_handle>     AssocEdgeMap;

  typedef std::pair<ISegment,ISegment>                   is_pair;

  struct lt_ssi_pair {

    public:
    lt_ssi_pair() {}
    bool operator()(const is_pair& s0, const is_pair& s1) const {
      if(s0.second == s1.second)
	return s0.first < s1.first;
      return s0.second < s1.second;
    }
  };

  typedef std::map<is_pair, event_iterator, lt_ssi_pair> EventHash;

  typedef std::multimap<Point_2, ISegment, lt_pnts_xy>   SegQueue;
  typedef typename SegQueue::iterator                    seg_iterator;
  typedef typename SegQueue::value_type                  ps_pair;

  event_iterator    event;
  Point_2           p_sweep;
  EventQueue        XS;
  seg_pair          sl,sh;
  compare_segs_at_sweepline SLcmp;
  SweepStatus       YS;
  Y2X               y2x;
  X2Y               x2y;
  X2iso             x2iso;
  Y2Y               y2y;
  SegQueue          SQ;
  IList             Internal;
  AssocEdgeMap      Edge_of;
  EventHash         IEvent;

  stl_seg_overlay_traits(const INPUT& in, OUTPUT& G, const GEOMETRY& k) 
    : its(in.first), ite(in.second), GO(G), K(k), 
    XS(compare_pnts_xy(K)), SLcmp(p_sweep,&sl,&sh,K), YS(SLcmp), 
    y2x(XS.end()), x2y(YS.end()), x2iso(0), y2y(&sl),
    SQ(lt_pnts_xy(K)), Edge_of(0)
  {}

  std::string dump_structures()
  { 
    std::ostringstream out; 
    out << "EventQueue:\n";
    typename EventQueue::const_iterator sit1;
    for(sit1 = XS.begin(); sit1 != XS.end(); ++sit1) 
      out << "  " << *sit1 << std::endl;

    out << "SegQueue:\n";
    typename SegQueue::const_iterator sit2;
    for(sit2 = SQ.begin(); sit2 != SQ.end(); ++sit2) 
      out << "  " << sit2->first << " " << sit2->second
          << " " << sit2->first << std::endl;

    out << "SweepStatus:\n";
    typename SweepStatus::iterator sit3;
    for( sit3 = YS.begin(); *sit3 != &sh; ++sit3 ) {
      int b = orientation(sit3, p_sweep);
      if(*sit3 == &sl) out << " 1";
      else if(*sit3 == &sh) out <<"-1";
      else if(b >= 0) out << " " << b;
      else out << b;
      out << " " << *sit3 << ": ";
      //      if(y2x[*sit3] != XS.end())
      //      	out << *y2x[*sit3];
      //      out << " | ";
      out << (*sit3)->first; 
      if(y2y[*sit3] != &sl)
	out << " y2y: " << y2y[*sit3];
      out << std::endl;
    }
    return out.str();
  }

  bool check_bundle(ss_iterator pred,
		    ss_iterator succ) {
    CGAL_NEF_TRACEN("check bundle");

    check_invariants();

    CGAL_assertion(*pred == &sl ||
		   orientation(pred, p_sweep) > 0);
    ++pred;
    while(*pred != *succ) {
      CGAL_assertion(orientation(pred, p_sweep) == 0);
      ss_iterator next(pred);
      ++next;
      if(*pred != &sl &&
	 *next != &sh) {
	bool b1 = orientation(pred, K.source((*next)->first)) == 0;
	b1 &= orientation(pred, K.target((*next)->first)) == 0;
	b1 &= orientation(next, K.source((*pred)->first)) == 0;
	b1 &= orientation(next, K.target((*pred)->first)) == 0;
	CGAL_warning(b1 == (y2y[*pred] == *next));
      }
      ++pred;
    }
    CGAL_assertion(*succ == &sh ||
		   orientation(succ, p_sweep) < 0);

    return true;
  }

  Point_2 source(ISegment is) const
  { return K.source(is->first); }
  Point_2 target(ISegment is) const
  { return K.target(is->first); }

  ITERATOR original(ISegment s) const
  { return s->second; }

  int orientation(ss_iterator sit, const Point_2& p) const
  { return K.orientation((*sit)->first,p); }

  bool collinear(ss_iterator sit1, ss_iterator sit2) const
  { if( *sit1 == &sl || 
	*sit1 == &sh || 
	*sit2 == &sl || 
	*sit2 == &sh) return false;
    Point_2 ps = source(*sit2), pt = target(*sit2);
    return ( orientation(sit1,ps)==0 &&
             orientation(sit1,pt)==0 );
  }

  event_iterator insertXS(const Point_2& p) {
    event_iterator upper = XS.upper_bound(p);
    if(upper == XS.begin())
      return XS.insert_before(upper, p);
    event_iterator pred = upper; --pred;
    if(K.compare_xy(*pred, p) == CGAL::SMALLER)
      return XS.insert_before(upper, p);
    return pred;
  }

  void compute_intersection(ss_iterator sit0)
  {    
    // Given an item |sit0| in the Y-structure compute the point of 
    // intersection with its successor and (if existing) insert it into 
    // the event queue and do all necessary updates.
    ss_iterator sit1 = sit0; ++sit1;
    CGAL_NEF_TRACEN("compute_intersection "<< *sit0 <<" "<< *sit1);
    if ( *sit0 == &sl || *sit1 == &sh ) return;
    const Segment_2& s0 = (*sit0)->first;
    const Segment_2& s1 = (*sit1)->first;
    int or0 = K.orientation(s0,K.target(s1));
    int or1 = K.orientation(s1,K.target(s0));
    if ( or0 <= 0 && or1 >= 0  ) { 
      event_iterator it = 
	IEvent[std::make_pair(*sit0, *sit1)];
      if(it == event_iterator()) {
	Point_2 q = K.intersection(s0,s1);
	event_iterator  er = insertXS(q); // only done if none existed!!!
	x2y[&*er] = sit0;
	y2x[*sit0] = er;
	CGAL_assertion(sit0 != YS.end());
      } else {
	CGAL_NEF_TRACEN("  intersection has been found previously");
	y2x[*sit0] = it;
      }
    }
  } 

  void initialize_structures()
  {
    /* INITIALIZATION
       - insert all vertices into the x-structure
       - insert sentinels into y-structure
       - exploit the fact that insert operations into the x-structure
         leave previously inserted points unchanged to achieve that
         any pair of endpoints $p$ and $q$ with |p == q| are identical
    */
    CGAL_NEF_TRACEN("initialize_structures");

    ITERATOR it_s;
    for ( it_s=its; it_s != ite; ++it_s ) {
      const Segment_2& s = *it_s;
      event_iterator it1, it2, upper;

      if(XS.empty())
	it1 = XS.insert(K.source(s));
      else
	it1 = insertXS(K.source(s));
      it2 = insertXS(K.target(s));
      
      if (it1 == it2) {
        if ( x2iso[&*it1] == 0 ) x2iso[&*it1] = new IsoList;
        x2iso[&*it1]->push_front(it_s);
        continue;  // ignore zero-length segments regarding YS
      }
          
      Point_2 p = *it1;
      Point_2 q = *it2;
          
      Segment_2 s1; 
      if ( K.compare_xy(p,q) < 0 ) 
        s1 = K.construct_segment(p,q);
      else
        s1 = K.construct_segment(q,p);
          
      Internal.push_back(seg_pair(s1,it_s));
      SQ.insert(ps_pair(K.source(s1),&Internal.back()));
    }

    // insert a lower and an upper sentinel segment to avoid special
    // cases when traversing the Y-structure
    YS.insert(&sl);
    YS.insert(&sh);
    CGAL_NEF_TRACEN("end of initialization\n");
  }


  bool event_exists() 
  { 
    if (!XS.empty()) { 
      // event is set at end of loop and in init
      event = XS.begin();
      p_sweep = *event;
      return true;
    }
    return false; 
  }

  void procede_to_next_event() 
  { XS.erase(event); }

  void process_event() 
  {
    CGAL_NEF_TRACEN("\n\n >>> process_event: "<<p_sweep);

    Vertex_handle v = GO.new_vertex(p_sweep);
    ss_iterator sit = x2y[&*event];
    ss_iterator sit_succ, sit_pred, sit_first, sit_pred_succ;

    if(sit == YS.end()) {
      CGAL_NEF_TRACEN("search for upper bound in YS");
      Segment_2 s_sweep = K.construct_segment(p_sweep,p_sweep);	
      seg_pair sp(s_sweep, ite);
      sit_succ = YS.upper_bound(&sp);
      sit = sit_succ;
      --sit;
      CGAL_NEF_TRACEN("upper bound: " << *sit_succ);
      if(*sit == &sl ||
	 orientation(sit, p_sweep) != 0) {
        sit_pred_succ = sit_succ;
	sit_pred = sit;
	sit = YS.end();
      }
    }
      
    /* |sit| is determined by upper bounding the search for the
       segment (p_sweep,p_sweep) and taking its predecessor.
       if the segment associated to |sit| contains |p_sweep| then
       there's a bundle of segments containing |p_sweep|.
       We compute the successor (|sit_succ)|) and 
       predecessor (|sit_pred|) items. */
    

    /* If no segments contain p_sweep then sit_pred and sit_succ are
       correctly set after the above locate operation, if a segment
       contains p_sweep sit_pred and sit_succ are set below when
       determining the bundle.*/
    
    if (sit != YS.end() ) { // sit->first is ending or passing segment
      CGAL_NEF_TRACEN("ending/passing segs " << *sit);
      sit_succ = sit; ++sit_succ;

      while ( y2x[*sit] == event ||
	      y2y[*sit] == *sit_succ ) { // overlapping
	sit = sit_succ;
	++sit_succ;
      }
      CGAL_NEF_TRACEN("ending/passing segs " << *sit);

      ss_iterator sit_last = sit;
      
      event_iterator xit = y2x[*sit_last];
      if (xit != XS.end() &&
	  *sit_last != &sl &&
	  *sit_succ != &sh) { 
	IEvent[std::make_pair(*sit_last, *sit_succ)] = xit;
      } 

      bool overlapping;
      do {
	ISegment s = *sit;
	ss_iterator sit_next(sit); --sit_next;

	if(*sit_next == &sl)
	  overlapping = false;
	else
	  overlapping = y2y[*sit_next] == *sit;
	Halfedge_handle e = Edge_of[*sit];
	if ( overlapping ) {
	  CGAL_NEF_TRACEN("overlapping segment "<<s);
	} else {
	  CGAL_NEF_TRACEN("connecting edge to node "<<s);
	  GO.link_as_target_and_append(v,e);
	  /* in this case we close the output edge |e| associated to
	     |sit| by linking |v| as its target and by appending the
	     twin edge to |v|'s adjacency list. */
	}
	GO.supporting_segment(e,original(s));
	if ( target(s) == p_sweep ) {
	  CGAL_NEF_TRACEN("ending segment "<<s);
	  if(overlapping)
	    y2y[*sit_next] = y2y[*sit];
	  YS.erase(sit);
	  GO.ending_segment(v,original(s));
	} else { // passing segment, take care of the node here!
	  CGAL_NEF_TRACEN("passing segment "<<s);
	  ss_iterator sst(sit); ++sst;
	  if(y2y[*sit] != *sst)
	    y2y[*sit] = &sl;
	  y2x[*sit] = XS.end();
	  GO.passing_segment(v,original(s));
	}

	ss_iterator sss(sit_next); ++sss;
	if(*sit_next != &sl)
	  overlapping |= (y2y[*sit_next] == *sss);

	sit = sit_next;
      }
      while (*sit != &sl && 
	     (y2x[*sit] == event || overlapping));
	     
      
      sit_pred = sit;
      sit_first = sit_pred;
      ++sit_first; // first item of the bundle
      sit_pred_succ = sit_first;

      CGAL_NEF_TRACE("event bundles between\n   "<<  *sit_succ);
      CGAL_NEF_TRACEN("\n   "<< *sit_pred);

      CGAL_assertion(check_bundle(sit_pred, sit_succ));

      while( *sit != *sit_succ) {
	ss_iterator sub_first(sit), 
	  sub_last(sub_first),
	  succ_sub_last(sub_last);
	++succ_sub_last;
	
	while(y2y[*sub_last] == *succ_sub_last) {
	  ++sub_last; ++succ_sub_last;
	}
	
	sit = sub_first;
	while(*sub_first != *sub_last) {
	  YS.swap(sub_first, sub_last);
	  std::swap(sub_first, sub_last);
	  ++sub_first;
	  if(*sub_first == *sub_last) break;
	  --sub_last;
	}
	++sit;
      }

      if(*sit_first != *sit_succ) {	
	ss_iterator sfirst(sit_pred); ++sfirst;
	ss_iterator slast(sit_succ); --slast;
	while(*sfirst != *slast) {
	  YS.swap(sfirst, slast);
	  std::swap(sfirst, slast);
	  ++sfirst;
	  if(*sfirst == *slast) break;
	  --slast;
	}
      }
    } // if (sit != ss_iterator() )
    
    //    CGAL_assertion( sit_pred != YS.end() );
    GO.halfedge_below(v,Edge_of[*sit_pred]);
    if ( x2iso[&*event] != 0 ) {
      IsoList* IL = x2iso[&*event];
      typename IsoList::const_iterator iso_it;
      for (iso_it = IL->begin(); iso_it != IL->end(); ++iso_it) 
	GO.trivial_segment(v,*iso_it);
      delete IL;
      x2iso[&*event] = 0;
    }
    
    ISegment next_seg = ISegment(); // to avoid /W4 warning
    seg_iterator next_it = SQ.begin();
    while ( next_it != SQ.end() && 
	    ( next_seg = next_it->second, p_sweep == source(next_seg)) ) {
      CGAL_NEF_TRACEN("inserting "<<next_seg);

      ss_iterator s_sit = YS.upper_bound(next_seg);
      ss_iterator p_sit(s_sit); --p_sit;
      
      sit = YS.insert_before(s_sit, next_seg);
      
      ss_iterator ys2_tmp = sit;
      ++ys2_tmp;
      CGAL_assertion(*s_sit == *ys2_tmp);

      if(*s_sit != &sh &&
	 orientation(s_sit, source(next_seg)) == 0 &&
	 orientation(s_sit, target(next_seg)) == 0) {
	y2y[*sit] = *s_sit;
      }

      if ( &sl != *p_sit &&
	   orientation(p_sit, source(next_seg) ) == 0 &&
	   orientation(p_sit, target(next_seg) ) == 0 ) {
	y2y[*p_sit] = *sit;
      }

      x2y[&*XS.find(target(next_seg))] = sit;
      GO.starting_segment(v,original(next_seg));
      
      // delete minimum and assign new minimum to next_seg
      SQ.erase(SQ.begin());
      next_it = SQ.begin();
    }
    
    // we insert new edge stubs, non-linked at target
    ss_iterator sit_curr = sit_succ, sit_prev = sit_succ;
    for( --sit_curr; *sit_curr != *sit_pred; 
	 sit_prev = sit_curr, --sit_curr ) {
      CGAL_NEF_TRACEN("checking outedge "<< *sit_curr <<"\n   "<< *sit_prev);
      if (y2y[*sit_curr] == *sit_prev) { // overlapping
	CGAL_assertion(collinear(sit_curr, sit_prev));
	Edge_of[*sit_curr] = Edge_of[*sit_prev];
      } else {
	CGAL_NEF_TRACEN("creating new edge ");
	Edge_of[*sit_curr] = GO.new_halfedge_pair_at_source(v);
      }
    }
    sit_first = sit_prev;
 
   event_iterator xxit = y2x[*sit_pred];
    if(xxit != XS.end() &&
       *sit_pred != &sl &&
       *sit_pred_succ != &sh) {
      IEvent[std::make_pair(*sit_pred, *sit_pred_succ)] = xxit;
      y2y[*sit_pred] = &sl;
      y2x[*sit_pred] = XS.end();
    } 
    
    CGAL_NEF_TRACEN("pred,succ = "<< *sit_pred <<" "<< *sit_succ);
    compute_intersection(sit_pred); 
    sit = sit_succ; --sit;
    if (*sit != *sit_pred)
      compute_intersection(sit); 
  }

  void complete_structures() {}
  void check_invariants() {CGAL_NEF_TRACEN("check_invariants\n"<<dump_structures());}
  void check_final() {}


}; // stl_seg_overlay_traits

} // namespace CGAL

#endif // !defined(CGAL_USE_LEDA_LIBRARY)

namespace CGAL {
#if defined(CGAL_USE_LEDA_LIBRARY)
#define Segment_overlay_traits leda_seg_overlay_traits
static const char* const sweepversion = "LEDA segment overlay sweep";
#else
#define Segment_overlay_traits stl_seg_overlay_traits
static const char* const sweepversion = "STL segment overlay sweep";
#endif
} // namespace CGAL

#include <CGAL/generic_sweep.h>
#endif // CGAL_SEGMENT_OVERLAY_TRAITS_H
