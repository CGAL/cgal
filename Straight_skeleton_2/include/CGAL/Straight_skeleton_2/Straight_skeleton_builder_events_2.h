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

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>

CGAL_BEGIN_NAMESPACE

namespace CGAL_SS_i
{

template<class SSkel_, class Traits_>
class Event_2 : public Ref_counted_base
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;
  
public:
 
  typedef Event_2<SSkel,Traits> Self ;

  typedef boost::intrusive_ptr<Self> SelfPtr ;

  typedef typename Traits::Point_2             Point_2 ;
  typedef typename Traits::FT                  FT ;
  typedef typename Traits::Seeded_trisegment_2 Seeded_trisegment_2 ;
  
  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;

  typedef Triedge<Halfedge_handle> Triedge ;
  
  enum Type { cEdgeEvent, cSplitEvent, cPseudoSplitEvent } ;

public:

  Event_2 ( Triedge const& aTriedge, Seeded_trisegment_2 const& aSTrisegment )
    :
     mTriedge    (aTriedge)
    ,mSTrisegment(aSTrisegment)
  {}

  virtual ~ Event_2() {}

  virtual Type type() const = 0 ;

  virtual Vertex_handle seed0() const = 0 ;
  virtual Vertex_handle seed1() const = 0 ;

  Triedge const&             triedge    () const { return mTriedge    ; }
  Seeded_trisegment_2 const& strisegment() const { return mSTrisegment; }
  Point_2 const&             point      () const { return mP          ; }
  FT                         time       () const { return mTime       ; }

  void SetTimeAndPoint( FT aTime, Point_2 const& aP ) { mTime = aTime ; mP = aP ; }

  friend std::ostream& operator<< ( std::ostream& ss, Self const& e )
  {
    ss << "[" ;
    e.dump(ss);
    ss << " p=(" << e.point().x() << "," << e.point().y() << ") t=" << e.time() << "] " 
       << trisegment_collinearity_to_string(e.strisegment().event().collinearity()) ;
    return ss ;
  }

protected :

  virtual void dump ( std::ostream& ss ) const
  {
    ss << mTriedge ;
  } ;

private :

  Triedge             mTriedge ;
  Seeded_trisegment_2 mSTrisegment ;
  Point_2             mP ;
  FT                  mTime ;
} ;

template<class SSkel_, class Traits_>
class Edge_event_2 : public Event_2<SSkel_,Traits_>
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;

  typedef Event_2<SSkel,Traits> Base ;

  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;

  typedef typename Base::Type                Type ;
  typedef typename Base::Triedge             Triedge ;
  typedef typename Base::Seeded_trisegment_2 Seeded_trisegment_2 ;

public:

  Edge_event_2 ( Triedge const&             aTriedge
               , Seeded_trisegment_2 const& aSTrisegment 
               , Vertex_handle              aLSeed
               , Vertex_handle              aRSeed
               )
    :
      Base(aTriedge,aSTrisegment)
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
  
  typedef typename Base::Type                Type ;
  typedef typename Base::Triedge             Triedge ;
  typedef typename Base::Seeded_trisegment_2 Seeded_trisegment_2 ;

public:

  Split_event_2 ( Triedge const&             aTriedge
                , Seeded_trisegment_2 const& aSTrisegment 
                , Vertex_handle              aSeed
                )
    :
      Base(aTriedge,aSTrisegment)
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

  typedef typename Base::Type                Type ;
  typedef typename Base::Triedge             Triedge ;
  typedef typename Base::Seeded_trisegment_2 Seeded_trisegment_2 ;

public:

  Pseudo_split_event_2 ( Triedge const&             aTriedge
                       , Seeded_trisegment_2 const& aSTrisegment 
                       , Vertex_handle              aSeed
                       , Vertex_handle              aOppositeNode
                       )
    :
      Base(aTriedge,aSTrisegment)
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
