// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_EVENTS_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_EVENTS_2_H 1

#include<ostream>

CGAL_BEGIN_NAMESPACE

namespace CGAL_SS_i
{

//
// This record encapsulates the defining contour halfedges for a node (both contour and skeleton)
//
template<class SSkel_>
struct Triedge
{
  typedef SSkel_   SSkel ;
  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  
  typedef Triedge<SSkel> Self ;
  
  Triedge() {}

  // Contour nodes (input polygon vertices) have just 2 defining contour edges    
  Triedge ( Halfedge_handle aE0, Halfedge_handle aE1 )
  {
    mE[0] = aE0 ;
    mE[1] = aE1 ;
  }              
  
  // Skeleton nodes (offset polygon vertices) have 3 defining contour edges    
  Triedge ( Halfedge_handle aE0, Halfedge_handle aE1 , Halfedge_handle aE2 )
  {
    mE[0] = aE0 ;
    mE[1] = aE1 ;
    mE[2] = aE2 ;
  }              
  
  Halfedge_handle e( unsigned idx ) const { CGAL_assertion(idx<3); return mE[idx]; }
  
  Halfedge_handle e0() const { return e(0); }
  Halfedge_handle e1() const { return e(1); }
  Halfedge_handle e2() const { return e(2); }
  
  bool is_valid() const { return e0() != e1() && e1() != e2() ; }
  
  // returns 1 if aE is one of the halfedges stored in this triedge, 0 otherwise.
  int contains ( Halfedge_handle aE ) const
  {
    return aE == e0() || aE == e1() || aE == e2() ? 1 : 0 ;
  }

  // Returns the number of common halfedges in the two triedges x and y
  static int CountInCommon( Self const& x, Self const& y )
  {
    return x.contains(y.e0()) + x.contains(y.e1()) + x.contains(y.e2()) ; 
  }
  
  // Returns true if the triedges store the same 3 halfedges (in any order)
  friend bool operator == ( Self const& x, Self const& y ) { return CountInCommon(x,y) == 3 ; }
  
  friend bool operator != ( Self const& x, Self const& y ) { return !(x==y) ; }
  
  friend Self operator & ( Self const& x, Self const& y )
  {
    return Self(x.e0(), x.e1(), ( x.e0() == y.e0() || x.e1() == y.e0() ) ? y.e1() : y.e0() ) ;
  }
  
  friend std::ostream& operator<< ( std::ostream& ss, Self const& t )
  {
    return ss << "{E" << t.e0()->id() << ",E" << t.e1()->id() << ",E" << t.e2()->id() << "}" ;
  }
  
  Halfedge_handle mE[3];
} ;

template<class SSkel_, class Traits_>
class Event_2 : public Ref_counted_base
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;
  
public:
 
  typedef Event_2<SSkel,Traits> Self ;

  typedef boost::intrusive_ptr<Self> SelfPtr ;

  typedef typename Traits::Point_2      Point_2 ;
  typedef typename Traits::FT           FT ;
  typedef typename Traits::Trisegment_2 Trisegment_2 ;

  typedef Triedge<SSkel> Triedge ;
  
  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;

  
  enum Type { cEdgeEvent, cSplitEvent, cPseudoSplitEvent } ;

public:

  Event_2 ( Triedge const& aTriedge, Trisegment_2 const& aTrisegment )
    :
     mTriedge   (aTriedge)
    ,mTrisegment(aTrisegment)
  {}

  virtual ~ Event_2() {}

  virtual Type type() const = 0 ;

  virtual Vertex_handle seed0() const = 0 ;
  virtual Vertex_handle seed1() const = 0 ;

  Triedge const&      triedge   () const { return mTriedge    ; }
  Trisegment_2 const& trisegment() const { return mTrisegment; }
  Point_2 const&      point     () const { return mP         ; }
  FT                  time      () const { return mTime      ; }

  void SetTimeAndPoint( FT aTime, Point_2 const& aP ) { mTime = aTime ; mP = aP ; }

  friend std::ostream& operator<< ( std::ostream& ss, Self const& e )
  {
    ss << "[" ;
    e.dump(ss);
    ss << " p=(" << e.point().x() << "," << e.point().y() << ") t=" << e.time() << "] " 
       << trisegment_collinearity_to_string(e.trisegment().collinearity()) ;
    return ss ;
  }

protected :

  virtual void dump ( std::ostream& ss ) const
  {
    ss << mTriedge ;
  } ;

private :

  Triedge      mTriedge ;
  Trisegment_2 mTrisegment ;
  Point_2      mP ;
  FT           mTime ;
} ;

template<class SSkel_, class Traits_>
class Edge_event_2 : public Event_2<SSkel_,Traits_>
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;

  typedef Event_2<SSkel,Traits> Base ;

  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;

  typedef typename Base::Type         Type ;
  typedef typename Base::Triedge      Triedge ;
  typedef typename Base::Trisegment_2 Trisegment_2 ;

public:

  Edge_event_2 ( Triedge const&      aTriedge
               , Trisegment_2 const& aTrisegment 
               , Vertex_handle       aLSeed
               , Vertex_handle       aRSeed
               )
    :
      Base(aTriedge,aTrisegment)
    , mLSeed(aLSeed)
    , mRSeed(aRSeed)
  {}

  virtual Type type() const { return this->cEdgeEvent ; }

  virtual Vertex_handle seed0() const { return mLSeed ; }
  virtual Vertex_handle seed1() const { return mRSeed ; }

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

template<class SSkel_, class Traits_>
class Split_event_2 : public Event_2<SSkel_,Traits_>
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;

  typedef Event_2<SSkel,Traits> Base ;

  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;
  
  typedef typename Base::Type         Type ;
  typedef typename Base::Triedge      Triedge ;
  typedef typename Base::Trisegment_2 Trisegment_2 ;

public:

  Split_event_2 ( Triedge const&      aTriedge
                , Trisegment_2 const& aTrisegment 
                , Vertex_handle       aSeed
                )
    :
      Base(aTriedge,aTrisegment)
    , mSeed(aSeed)
  {}

  virtual Type type() const { return this->cSplitEvent ; }

  virtual Vertex_handle seed0() const { return mSeed ; }
  virtual Vertex_handle seed1() const { return mSeed ; }

  void set_opposite_rnode( Vertex_handle aOppR ) { mOppR = aOppR ; }
  
  Vertex_handle opposite_rnode() const { return mOppR ; }
  
private :

  virtual void dump ( std::ostream& ss ) const
  {
    this->Base::dump(ss);
    ss << " (Seed=" << mSeed->id() << " OppBorder=" << this->triedge().e2()->id() << ')' ;
  }

private :

  Vertex_handle mSeed ;
  Vertex_handle mOppR ;
} ;

template<class SSkel_, class Traits_>
class Pseudo_split_event_2 : public Event_2<SSkel_,Traits_>
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;
  
  typedef Event_2<SSkel,Traits> Base ;

  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;

  typedef typename Base::Type         Type ;
  typedef typename Base::Triedge      Triedge ;
  typedef typename Base::Trisegment_2 Trisegment_2 ;

public:

  Pseudo_split_event_2 ( Triedge const&      aTriedge
                       , Trisegment_2 const& aTrisegment 
                       , Vertex_handle       aSeed
                       , Vertex_handle       aOppositeNode
                       )
    :
      Base(aTriedge,aTrisegment)
    , mSeed(aSeed)
    , mOppNode(aOppositeNode)
  {}

  virtual Type type() const { return this->cPseudoSplitEvent ; }

  virtual Vertex_handle seed0() const { return mSeed ; }
  virtual Vertex_handle seed1() const { return mOppNode ; }

private :

  virtual void dump ( std::ostream& ss ) const
  {
    this->Base::dump(ss);
    ss << " (Seed=" << mSeed->id() << " OppNode=" << mOppNode->id() << ')' ;
  }

private :

  Vertex_handle mSeed   ;
  Vertex_handle mOppNode ;
} ;

}

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_EVENTS_2_H //
// EOF //
