// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
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
// file          : include/CGAL/straight_skeleton_builder_events_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_EVENTS_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_EVENTS_2_H 1

#include<ostream>

CGAL_BEGIN_NAMESPACE

template<class R>
class Straight_skeleton_builder_event_2 : public Ref_counted_base
{
public:
  
  typedef typename R::Rep Rep ;
   
  typedef typename Rep::Point_2 Point_2 ;
  typedef typename Rep::FT      FT ;
  
  typedef typename R::Halfedge_handle Halfedge_handle ;
  
  enum Type { cEdgeEvent, cSplitEvent } ;
  
public:

  Straight_skeleton_builder_event_2 (  Halfedge_handle aBorderA
                                     , Halfedge_handle aBorderB
                                     , Halfedge_handle aBorderC
                                    )
    :
     mBorderA(aBorderA)
    ,mBorderB(aBorderB)
    ,mBorderC(aBorderC)
  {}  
  
  virtual ~ Straight_skeleton_builder_event_2() {}
  
  virtual Type type() const = 0 ;
  
  Halfedge_handle border_a() const { return mBorderA ; }
  Halfedge_handle border_b() const { return mBorderB ; }
  Halfedge_handle border_c() const { return mBorderC ; }
  Point_2         point   () const { return mP ; }
  FT              time    () const { return mTime ; }
  
  void SetTimeAndPoint( FT aTime, Point_2 const& aP ) { mTime = aTime ; mP = aP ; }
  
  friend std::ostream& operator<< ( std::ostream& ss
                                   ,Straight_skeleton_builder_event_2<R> const& e 
                                  )
  {
    ss << "[" ;
    e.dump(ss);
    ss << " p=(" << e.point().x() << "," << e.point().y() << ") t=" << e.time() << "]" ;
    return ss ;
  }
     
protected :

  virtual void dump ( std::ostream& ss ) const 
  {
    ss << "{E" << mBorderA->id() << ",E" << mBorderB->id() << ",E" << mBorderC->id() << '}' ;   
  } ;
  
private :

  Halfedge_handle mBorderA ;
  Halfedge_handle mBorderB ;
  Halfedge_handle mBorderC ;
  Point_2         mP ;
  FT              mTime ;
} ;

template<class R>
class Straight_skeleton_builder_edge_event_2 : public Straight_skeleton_builder_event_2<R>
{
  
  typedef Straight_skeleton_builder_event_2<R> Base ;
    
  typedef typename R::Rep Rep ;
  
  typedef typename Rep::Point_2 Point_2 ;
  typedef typename Rep::FT      FT ;
  
  typedef typename R::Halfedge_handle Halfedge_handle ;
  typedef typename R::Vertex_handle   Vertex_handle ;
  
  typedef typename Base::Type Type ;
  
public:

  Straight_skeleton_builder_edge_event_2 (  Halfedge_handle aBorderA
                                          , Halfedge_handle aBorderB
                                          , Halfedge_handle aBorderC
                                          , Vertex_handle   aLSeed
                                          , Vertex_handle   aRSeed
                                          )
    :
      Base(aBorderA,aBorderB,aBorderC)
    , mLSeed(aLSeed)
    , mRSeed(aRSeed)
  {}  
  
  virtual Type type() const { return cEdgeEvent ; }
  
  Vertex_handle left_seed () const { return mLSeed ; }
  Vertex_handle right_seed() const { return mRSeed ; }
    
private :

  virtual void dump ( std::ostream& ss ) const 
  {
    this->Base::dump(ss);
    ss << " (LSeed=" << mLSeed->id() << " RSeed=" << mRSeed->id() << ')' ;
  }
  
private :

  Vertex_handle mLSeed ;
  Vertex_handle mRSeed ;
} ;

template<class R>
class Straight_skeleton_builder_split_event_2 : public Straight_skeleton_builder_event_2<R>
{
  
  typedef Straight_skeleton_builder_event_2<R> Base ;
  
  typedef Straight_skeleton_builder_split_event_2<R> Self ;
    
  typedef typename R::Rep Rep ;
  
  typedef typename Rep::Point_2 Point_2 ;
  typedef typename Rep::FT      FT ;
  
  typedef typename R::Halfedge_handle Halfedge_handle ;
  typedef typename R::Vertex_handle   Vertex_handle ;
  typedef typename Base::Type Type ;  
public:

  Straight_skeleton_builder_split_event_2 (  Halfedge_handle aBorderA
                                           , Halfedge_handle aBorderB
                                           , Halfedge_handle aBorderC
                                           , Vertex_handle   aSeed
                                           , Halfedge_handle aOppositeBorder
                                         )
    :
      Base(aBorderA,aBorderB,aBorderC)
    , mSeed(aSeed)
    , mOppositeBorder(aOppositeBorder)
    , mNext(0)
  {}  
  
  virtual Type type() const { return cSplitEvent ; }
  
  Vertex_handle seed() const { return mSeed ; }
  
  Halfedge_handle opposite_border() const { return mOppositeBorder ; }

  //
  // A "Vertex" event is a set of simultaneous split events.
  // It is represented as a linked list of split events.
  //
  Self* next() { return mNext ; }
  void set_next ( Self* aNext ) { mNext = aNext ; }
      
private :

  virtual void dump ( std::ostream& ss ) const 
  {
    this->Base::dump(ss);
    ss << " (Seed=" << mSeed->id() << " OppBorder=" << mOppositeBorder->id() << ')' ;
  }
  
private :

  Vertex_handle   mSeed ;
  Halfedge_handle mOppositeBorder ;
  Self*           mNext ;
} ;

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_EVENTS_2_H //
// EOF //
