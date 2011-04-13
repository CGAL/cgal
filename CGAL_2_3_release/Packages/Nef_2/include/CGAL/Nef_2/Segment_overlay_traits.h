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
// file          : include/CGAL/Nef_2/Segment_overlay_traits.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// source        : nef_2d/Segment_overlay.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: generic segment intersection sweep
// ============================================================================
#ifndef CGAL_SEGMENT_OVERLAY_TRAITS_H
#define CGAL_SEGMENT_OVERLAY_TRAITS_H

#include <assert.h>
#undef _DEBUG
#define _DEBUG 23
#include <CGAL/Nef_2/debug.h>

//#define INCLUDEBOTH
#if defined(CGAL_USE_LEDA) || defined(INCLUDEBOTH)
#include <LEDA/tuple.h>
#include <LEDA/slist.h>
#include <LEDA/list.h>
#include <LEDA/map.h>
#include <LEDA/map2.h>
#include <LEDA/sortseq.h>
#include <LEDA/p_queue.h>
#include <utility>
#include <strstream>

namespace CGAL {
#ifdef _DEBUG
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

  class cmp_segs_at_sweepline : public leda_cmp_base<ISegment>
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
     else ASSERT(0,"compare error in sweep.");
     if ( s || K.is_degenerate(s1) || K.is_degenerate(s2) ) 
       return s;
    
     s = K.orientation(s2,K.target(s1));
     if (s==0) return ( is1 - is2 );
     // overlapping segments are not equal
     return s;
   }
  };

  struct cmp_pnts_xy : public leda_cmp_base<Point_2>
  { const GEOMETRY& K;
  public:
   cmp_pnts_xy(const GEOMETRY& k) : K(k) {}
   int operator()(const Point_2& p1, const Point_2& p2) const
   { return K.compare_xy(p1,p2); }
  };


  typedef leda_sortseq<Point_2,seq_item>      EventQueue; 
  typedef leda_sortseq<ISegment,seq_item>     SweepStatus;
  typedef leda_p_queue<Point_2,ISegment>      SegQueue; 
  typedef leda_map<seq_item,Halfedge_handle>  AssocEdgeMap;
  typedef leda_slist<ITERATOR>                IsoList;
  typedef leda_map<seq_item, IsoList* >       AssocIsoMap;
  typedef leda_map2<ISegment,ISegment,seq_item> EventHash;

    seq_item                   event;
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
    std::ostrstream out;
    out << "SQ= ";
    pq_item pqit;
    forall_items(pqit,SQ) {
      if (SQ.prio(pqit)==XS.key(XS.succ(XS.min_item()))) 
      { out << SQ.inf(pqit)->first(); }
      pqit = SQ.next_item(pqit);
    }
    seq_item sit;
    out << "\nXS=\n";
    forall_items(sit,XS)
      out << "  " << XS.key(sit) << " " << XS.inf(sit) 
          <<std::endl;
    out << "YS=\n";
    for( sit = YS.max_item(); sit; sit=YS.pred(sit) )
      out << "  "<<YS.key(sit)->first()<<" "<<YS.inf(sit)<<std::endl;
    out << '\0';
    leda_string res(out.str()); out.freeze(0); return res;
  }


  Point_2 source(ISegment is) const
  { return K.source(is->first()); }
  Point_2 target(ISegment is) const
  { return K.target(is->first()); }
  ITERATOR original(ISegment s) const
  { return s->second(); }

  int orientation(seq_item sit, const Point_2& p) const
  { return K.orientation(YS.key(sit)->first(),p); }

  bool collinear(seq_item sit1, seq_item sit2) const
  { Point_2 ps = source(YS.key(sit2)), pt = target(YS.key(sit2));
    return ( orientation(sit1,ps)==0 &&
             orientation(sit1,pt)==0 );
  }


  void compute_intersection(seq_item sit0)
  {    
    seq_item sit1 = YS.succ(sit0);
    if ( sit0 == YS.min_item() || sit1 == YS.max_item() ) return;
    ISegment s0 = YS.key(sit0);
    ISegment s1 = YS.key(sit1);
    int or0 = K.orientation(s0->first(),target(s1));
    int or1 = K.orientation(s1->first(),target(s0));
    if ( or0 <= 0 && or1 >= 0  ) { 
      seq_item it = IEvent(YS.key(sit0),YS.key(sit1));
      if ( it==0 ) {
        Point_2 q = K.intersection(s0->first(),s1->first());
        it = XS.insert(q,sit0);
      }
      YS.change_inf(sit0, it);
    }
  } 

  void initialize_structures()
  {
    TRACEN("initialize_structures");
    ITERATOR it_s;  
    for ( it_s=its; it_s != ite; ++it_s ) {
      Segment_2 s = *it_s;
      seq_item it1 = XS.insert( K.source(s), seq_item(nil));
      seq_item it2 = XS.insert( K.target(s), seq_item(nil));
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
    YS.insert(&sl,seq_item(nil));
    YS.insert(&sh,seq_item(nil));
    TRACEN("end of initialization\n"<<YS.size());
  }

  bool event_exists() 
  { 
    if (!XS.empty()) { 
      // event is set at end of loop and in init
      event = XS.min();
      p_sweep = XS.key(event);
      return true;
    }
    return false; 
  }

  void procede_to_next_event() 
  { XS.del_item(event); }

  void process_event() 
  {
    TRACEN("\n\n >>> process_event: "<<p_sweep<<" "<<XS[event]<<" "<<event);

    Vertex_handle v = GO.new_vertex(p_sweep);
    seq_item sit = XS.inf(event);
        
      seq_item sit_succ(0), sit_pred(0), sit_pred_succ(0), sit_first(0);
      if (sit == nil) 
        {
          Segment_2 s_sweep = K.construct_segment(p_sweep,p_sweep);
          seg_pair sp(s_sweep,ITERATOR());
          sit_succ = YS.locate( &sp );
          if ( sit_succ != YS.max_item() && 
               orientation(sit_succ,p_sweep) == 0 ) 
            sit = sit_succ;
          else  {
            sit_pred = YS.pred(sit_succ);
            sit_pred_succ = sit_succ;
          }
          TRACEN("looked up p_sweep "<<PIS(YS.key(sit_succ)));
        }



      /* If no segment contains p_sweep then sit_pred and sit_succ are
         correctly set after the above locate operation, if a segment
         contains p_sweep sit_pred and sit_succ are set below when
         determining the bundle.*/

      if (sit != nil) { // key(sit) is an ending or passing segment
        TRACEN("ending/passing segs");
        while ( YS.inf(sit) == event ||
                YS.inf(sit) == YS.succ(sit) ) // overlapping
          sit = YS.succ(sit);
        sit_succ = YS.succ(sit); 
        seq_item sit_last = sit;

        seq_item xit = YS.inf(sit_last);
        if (xit) { 
          ISegment s1 = YS.key(sit_last);
          ISegment s2 = YS.key(sit_succ);
          IEvent(s1,s2) = xit;
            TRACEN("hashing "<<PIS(s1)<<PIS(s2)<<xit);
        } 
          
        bool overlapping;
        do {
          ISegment s = YS.key(sit);
          seq_item sit_next = YS.pred(sit);
          overlapping = (YS.inf(sit_next) == sit);
          Halfedge_handle e = Edge_of[sit];
          if ( !overlapping ) {
            TRACEN("connecting edge to node "<<PIS(s)<<" "<<sit);
            GO.link_as_target_and_append(v,e);
          }
          GO.supporting_segment(e,original(s));
          if ( target(s) == p_sweep ) { // ending segment
              TRACEN("ending segment "<<PIS(s));
            if ( overlapping ) YS.change_inf(sit_next,YS.inf(sit));
            YS.del_item(sit);
            GO.ending_segment(v,original(s));
          } else {  // passing segment
              TRACEN("passing segment "<<PIS(s));
            if ( YS.inf(sit) != YS.succ(sit) ) 
              YS.change_inf(sit, seq_item(0));
            GO.passing_segment(v,original(s));
          }
          sit = sit_next;
        } 
        while ( YS.inf(sit) == event || overlapping ||
                YS.inf(sit) == YS.succ(sit) );
                  
        sit_pred = sit;
        sit_first = sit_pred_succ = YS.succ(sit_pred); // first item of bundle

        TRACE("event bundles between\n   "<<PIS(YS.key(sit_succ)));
        TRACEN("\n   "<<PIS(YS.key(sit_pred)));

        while ( sit != sit_succ ) {
          seq_item sub_first = sit;
          seq_item sub_last  = sub_first;
                            
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

    assert(sit_pred);
    GO.halfedge_below(v,Edge_of[sit_pred]);
    if ( Isos_of[event] != 0 ) {
      const IsoList& IL = *(Isos_of[event]);
      slist_item iso_it;
      for (iso_it = IL.first(); iso_it; iso_it=IL.succ(iso_it) ) 
        GO.trivial_segment(v,IL[iso_it] );
      delete (Isos_of[event]); // clean up the list
    }


    ISegment next_seg;
    pq_item next_it = SQ.find_min();
    while ( next_it && 
            (next_seg = SQ.inf(next_it), p_sweep == source(next_seg)) ) {
      seq_item s_sit = YS.locate_succ(next_seg);
      seq_item p_sit = YS.pred(s_sit);

      TRACEN("inserting "<<PIS(next_seg)<<" at "<<PIS(YS.key(s_sit))); 
      if ( YS.max_item() != s_sit &&
           orientation(s_sit, source(next_seg) ) == 0 &&
           orientation(s_sit, target(next_seg) ) == 0 )
        sit = YS.insert_at(s_sit, next_seg, s_sit);
      else 
        sit = YS.insert_at(s_sit, next_seg, seq_item(nil));
      assert(YS.succ(sit)==s_sit);

      if ( YS.min_item() != p_sit &&
           orientation(p_sit, source(next_seg) ) == 0 &&
           orientation(p_sit, target(next_seg) ) == 0 )
        YS.change_inf(p_sit, sit);
      assert(YS.succ(p_sit)==sit);
                 
      XS.insert(target(next_seg), sit);
      GO.starting_segment(v,original(next_seg));
                 
      // delete minimum and assign new minimum to next_seg
      SQ.del_min();    
      next_it = SQ.find_min();
    }

    for( seq_item sitl = YS.pred(sit_succ); sitl != sit_pred; 
         sitl = YS.pred(sitl) ) {
      if ( YS.inf(sitl) != YS.succ(sitl) ) { // non-overlapping
        TRACEN("non-overlapping "<<PIS(YS.key(sitl))<<" "<<sitl);
        Edge_of[sitl] = GO.new_halfedge_pair_at_source(v);
      } else {
        TRACEN("overlapping "<<PIS(YS.key(sitl)));
        Edge_of[sitl] = Edge_of[ YS.succ(sitl) ];
      }
    }
    sit_first = YS.succ(sit_pred);


    assert(sit_pred); assert(sit_pred_succ);
    seq_item xit = YS.inf(sit_pred);
    if ( xit ) { 
      ISegment s1 = YS.key(sit_pred);
      ISegment s2 = YS.key(sit_pred_succ);
      IEvent(s1,s2) = xit;
        TRACEN("hashing "<<PIS(s1)<<PIS(s2)<<xit);
      YS.change_inf(sit_pred, seq_item(0));
    }
          
    compute_intersection(sit_pred); 
    sit = YS.pred(sit_succ);
    if (sit != sit_pred)
      compute_intersection(sit);


  }

  void complete_structures() {}
  void check_invariants() {TRACEN("check_invariants\n"<<dump_structures());}
  void check_final() {}

}; // leda_seg_overlay_traits

} // namespace CGAL

#endif // defined(CGAL_USE_LEDA) || defined(INCLUDEBOTH)
#if !defined(CGAL_USE_LEDA) || defined(INCLUDEBOTH)
#include <list>
#include <map>
#include <string>
#include <strstream>

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
  { const Point_2& p;
    ISegment s_bottom, s_top; // sentinel segments
    const GEOMETRY& K;
  public:
   lt_segs_at_sweepline(const Point_2& pi, 
     ISegment s1, ISegment s2, const GEOMETRY& k) : 
     p(pi), s_bottom(s1), s_top(s2), K(k) {}
   lt_segs_at_sweepline(const lt_segs_at_sweepline& lt) :
     p(lt.p), s_bottom(lt.s_bottom), s_top(lt.s_top), K(lt.K) {}

   bool operator()(const ISegment& is1, const ISegment& is2) const
   { 
     if ( is2 == s_top || is1 == s_bottom ) return true;
     if ( is1 == s_top || is2 == s_bottom ) return false;
     if ( is1 == is2 ) return false;
     // Precondition: p is contained in s1 or s2. 
     const Segment_2& s1 = is1->first;
     const Segment_2& s2 = is2->first;
     int s = 0;
     if ( K.orientation(s1,p) == 0 )      
       s =   K.orientation(s2,p);
     else if ( K.orientation(s2,p) == 0 ) 
       s = - K.orientation(s1,p);
     else ASSERT(0,"compare error in sweep.");
     if ( s || K.is_degenerate(s1) || K.is_degenerate(s2) ) 
       return ( s < 0 );
    
     s = K.orientation(s2,K.target(s1));
     if (s==0) return ( is1 - is2 ) < 0;
     // overlapping segments are not equal
     return ( s < 0 );
   }
  };

  struct lt_pnts_xy 
  { const GEOMETRY& K;
  public:
   lt_pnts_xy(const GEOMETRY& k) : K(k) {}
   lt_pnts_xy(const lt_pnts_xy& lt) : K(lt.K) {}
   int operator()(const Point_2& p1, const Point_2& p2) const
   { return K.compare_xy(p1,p2) < 0; }
  };


  typedef std::map<ISegment, Halfedge_handle, lt_segs_at_sweepline> 
                                                         SweepStatus;
  typedef typename SweepStatus::iterator                 ss_iterator;
  typedef typename SweepStatus::value_type               ss_pair;

  typedef std::list<ITERATOR>                            IsoList;
  typedef std::map<Point_2, IsoList*, lt_pnts_xy>        EventQueue;
  typedef typename EventQueue::iterator                  event_iterator;
  typedef typename EventQueue::value_type                event_pair;

  typedef std::multimap<Point_2, ISegment, lt_pnts_xy>   SegQueue;
  typedef typename SegQueue::iterator                    seg_iterator;
  typedef typename SegQueue::value_type                  ps_pair;

  event_iterator    event;
  Point_2           p_sweep;
  EventQueue        XS;
  seg_pair          sl,sh;
  SweepStatus       YS;
  SegQueue          SQ;
  IList             Internal;

  stl_seg_overlay_traits(const INPUT& in, OUTPUT& G, 
    const GEOMETRY& k) : 
    its(in.first), ite(in.second), GO(G), K(k), 
    XS(lt_pnts_xy(K)), YS(lt_segs_at_sweepline(p_sweep,&sl,&sh,K)),
    SQ(lt_pnts_xy(K)) {}


  std::string dump_structures() const
  { 
    std::ostrstream out; 
    out << "EventQueue:\n";
    typename EventQueue::const_iterator sit1;
    for(sit1 = XS.begin(); sit1 != XS.end(); ++sit1) 
      out << "  " << sit1->first << std::endl;

    out << "SegQueue:\n";
    typename SegQueue::const_iterator sit2;
    for(sit2 = SQ.begin(); sit2 != SQ.end(); ++sit2) 
      out << "  " << sit2->first << " " << sit2->second
          << " " << sit2->first << std::endl;

    out << "SweepStatus:\n";
    typename SweepStatus::const_iterator sit3;
    for( sit3 = YS.begin(); sit3 != YS.end(); ++sit3 )
      out << sit3->first << " " << &*(sit3->second) << std::endl;
    std::string res(out.str()); out.freeze(0); return res;
  }

  Point_2 source(ISegment is) const
  { return K.source(is->first); }
  Point_2 target(ISegment is) const
  { return K.target(is->first); }

  ITERATOR original(ISegment s) const
  { return s->second; }

  int orientation(ss_iterator sit, const Point_2& p) const
  { return K.orientation(sit->first->first,p); }

  bool collinear(ss_iterator sit1, ss_iterator sit2) const
  { Point_2 ps = source(sit2->first), pt = target(sit2->first);
    return ( orientation(sit1,ps)==0 &&
             orientation(sit1,pt)==0 );
  }

  void compute_intersection(ss_iterator sit0)
  {    
    // Given an item |sit0| in the Y-structure compute the point of 
    // intersection with its successor and (if existing) insert it into 
    // the event queue and do all necessary updates.
    ss_iterator sit1 = sit0; ++sit1;
    TRACEN("compute_intersection "<<sit0->first<<" "<<sit1->first);
    if ( sit0 == YS.begin() || sit1 == --YS.end() ) return;
    const Segment_2& s0 = sit0->first->first;
    const Segment_2& s1 = sit1->first->first;
    int or0 = K.orientation(s0,K.target(s1));
    int or1 = K.orientation(s1,K.target(s0));
    if ( or0 <= 0 && or1 >= 0  ) { 
      Point_2 q = K.intersection(s0,s1);
      XS.insert(event_pair(q,0)); // only done if none existed!!!
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
    TRACEN("initialize_structures");
    
    ITERATOR it_s;  
    for ( it_s=its; it_s != ite; ++it_s ) {
      const Segment_2& s = *it_s;
      event_iterator it1 = (XS.insert(event_pair(K.source(s),0))).first;
      event_iterator it2 = (XS.insert(event_pair(K.target(s),0))).first;
      // note that the STL only inserts if key is not yet in XS
      if (it1 == it2) { 
        if ( it1->second == 0 ) it1->second = new IsoList;
        it1->second->push_front(it_s);
        continue;  // ignore zero-length segments regarding YS
      }
          
      Point_2 p = it1->first;
      Point_2 q = it2->first;
          
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
    YS.insert(ss_pair(&sl,Halfedge_handle()));
    YS.insert(ss_pair(&sh,Halfedge_handle()));
    TRACEN("end of initialization\n");
  }


  bool event_exists() 
  { 
    if (!XS.empty()) { 
      // event is set at end of loop and in init
      event = XS.begin();
      p_sweep = event->first;
      return true;
    }
    return false; 
  }

  void procede_to_next_event() 
  { XS.erase(event); }

  void process_event() 
  {
    TRACEN("\n\n >>> process_event: "<<p_sweep);

    Vertex_handle v = GO.new_vertex(p_sweep);
    ss_iterator sit_succ, sit_pred, sit_first, sit;
    Segment_2 s_sweep = K.construct_segment(p_sweep,p_sweep);
    seg_pair sp(s_sweep,ITERATOR());
    sit_succ = YS.upper_bound(&sp);
    sit = sit_succ; --sit;

    
      /* |sit| is determined by upper bounding the search for the
         segment (p_sweep,p_sweep) and taking its predecessor.
         if the segment associated to |sit| contains |p_sweep| then
         there's a bundle of segments containing |p_sweep|.
         We compute the successor (|sit_succ)|) and 
         predecessor (|sit_pred|) items. */
      
       if ( sit == YS.begin() || orientation(sit,p_sweep) != 0 ) {
        sit_pred = sit;
        sit = YS.end();
      }

      /* If no segments contain p_sweep then sit_pred and sit_succ are
         correctly set after the above locate operation, if a segment
         contains p_sweep sit_pred and sit_succ are set below when
         determining the bundle.*/

      if ( sit != YS.end() ) { // sit->first is ending or passing segment
        // Determine upper bundle item:
        TRACEN("ending/passing segs");
      
        /* Walk down until |sit_pred|, close edges for all segments 
           in the bundle, delete all segments in the bundle, but 
           reinsert the continuing ones */

        std::list<ISegment> L_tmp;
        bool overlapping;
        do {
          ISegment s = sit->first;
          ss_iterator sit_next(sit); --sit_next;
          overlapping = (sit_next != YS.begin()) && collinear(sit,sit_next);
          Halfedge_handle e = sit->second;
          if ( overlapping ) {
            TRACEN("overlapping segment "<<s);
          } else {
            TRACEN("connecting edge to node "<<s);
            GO.link_as_target_and_append(v,e);
            /* in this case we close the output edge |e| associated to 
               |sit| by linking |v| as its target and by appending the 
               twin edge to |v|'s adjacency list. */
          }
          GO.supporting_segment(e,original(s));

          if ( target(s) == p_sweep ) {
              TRACEN("ending segment "<<s);
            GO.ending_segment(v,original(s));
          } else { // passing segment, take care of the node here!
              TRACEN("passing segment "<<s);
            L_tmp.push_back(s);
            GO.passing_segment(v,original(s));
           }
          sit = sit_next;
        } 
        while ( sit != YS.begin() && orientation(sit,p_sweep) == 0 );
              
        sit_pred = sit_first = sit;
        ++sit_first; // first item of the bundle

        TRACE("event bundles between\n   "<<sit_succ->first);
        TRACEN("\n   "<<sit_pred->first);

        /* Interfaceproposition for next chunk:
           - succ(sit_pred) == sit_first == sit_succ
           - bundle not empty: sit_first != sit_succ
        */

        // delete and reinsert the continuing bundle
        YS.erase(sit_first,sit_succ);
        typename std::list<ISegment>::const_iterator lit;
        for ( lit = L_tmp.begin(); lit != L_tmp.end(); ++lit ) {
          YS.insert(sit_pred,ss_pair(*lit,Halfedge_handle()));
        }
      } // if (sit != ss_iterator() )


      assert( sit_pred != YS.end() );
      GO.halfedge_below(v,sit_pred->second);
      if ( event->second != 0 ) {
        const IsoList& IL = *(event->second);
        typename IsoList::const_iterator iso_it;
        for (iso_it = IL.begin(); iso_it != IL.end(); ++iso_it) 
          GO.trivial_segment(v,*iso_it);
        delete (event->second);
      }


      ISegment next_seg;
      seg_iterator next_it = SQ.begin();
      while ( next_it != SQ.end() && 
              ( next_seg = next_it->second, p_sweep == source(next_seg)) ) {
        TRACEN("inserting "<<next_seg);
        YS.insert(ss_pair(next_seg,Halfedge_handle()));
        GO.starting_segment(v,original(next_seg));
        // delete minimum and assign new minimum to next_seg
        SQ.erase(SQ.begin());
        next_it = SQ.begin();
      }
      // we insert new edge stubs, non-linked at target
      ss_iterator sit_curr = sit_succ, sit_prev = sit_succ;
      for( --sit_curr; sit_curr != sit_pred; 
           sit_prev = sit_curr, --sit_curr ) {
        TRACEN("checking outedge "<<sit_curr->first<<"\n   "<<sit_prev->first);
        if ( sit_curr != YS.begin() && sit_prev != --YS.end() &&
             collinear(sit_curr,sit_prev) ) // overlapping
          sit_curr->second = sit_prev->second;
        else {
          TRACEN("creating new edge");
          sit_curr->second = GO.new_halfedge_pair_at_source(v);
        }
      }
      sit_first = sit_prev;


      // compute possible intersections between |sit_pred| and its
      // successor and |sit_succ| and its predecessor
      TRACEN("pred,succ = "<<sit_pred->first<<" "<<sit_succ->first);
      compute_intersection(sit_pred); 
      sit = sit_succ; --sit;
      if (sit != sit_pred)
        compute_intersection(sit);

  }

  void complete_structures() {}
  void check_invariants() {TRACEN("check_invariants\n"<<dump_structures());}
  void check_final() {}


}; // stl_seg_overlay_traits

} // namespace CGAL

#endif // !defined(CGAL_USE_LEDA) || defined(INCLUDEBOTH)

namespace CGAL {
#ifdef CGAL_USE_LEDA
#define Segment_overlay_traits leda_seg_overlay_traits
static const char* sweepversion = "LEDA segment overlay sweep";
#else
#define Segment_overlay_traits stl_seg_overlay_traits
static const char* sweepversion = "STL segment overlay sweep";
#endif
} // namespace CGAL

#include <CGAL/generic_sweep.h>
#endif // CGAL_SEGMENT_OVERLAY_TRAITS_H

