// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_EVENTS_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_EVENTS_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include<ostream>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>

namespace CGAL {

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

  typedef typename Traits::Point_2          Point_2 ;
  typedef typename Traits::FT               FT ;
  typedef typename Traits::Trisegment_2_ptr Trisegment_2_ptr ;
  
  typedef typename SSkel::Vertex                Vertex ;
  typedef typename SSkel::Halfedge_handle       Halfedge_handle ;
  typedef typename SSkel::Halfedge_const_handle Halfedge_const_handle ;
  typedef typename SSkel::Vertex_handle         Vertex_handle ;
  
  typedef CGAL_SS_i::Triedge<Halfedge_handle> Triedge ;
  
  enum Type { cEdgeEvent, cSplitEvent, cPseudoSplitEvent } ;

public:

  Event_2 ( Triedge const& aTriedge, Trisegment_2_ptr const& aTrisegment )
    :
     mTriedge   (aTriedge)
    ,mTrisegment(aTrisegment)
  {}

  virtual ~ Event_2() {}

  virtual Type type() const = 0 ;

  virtual Vertex_handle seed0() const = 0 ;
  virtual Vertex_handle seed1() const = 0 ;

  Triedge const&          triedge   () const { return mTriedge   ; }
  Trisegment_2_ptr const& trisegment() const { return mTrisegment; }
  Point_2 const&          point     () const { return mP         ; }
  FT                      time      () const { return mTime      ; }

  void SetTimeAndPoint( FT aTime, Point_2 const& aP ) { mTime = aTime ; mP = aP ; }

  friend std::ostream& operator<< ( std::ostream& ss, Self const& e )
  {
    ss << "[" ;
    e.dump(ss);
    ss << " p=(" << e.point().x() << "," << e.point().y() << ") t=" << e.time() << "] " 
       << trisegment_collinearity_to_string(e.trisegment()->collinearity()) ;
    return ss ;
  }

protected :

  virtual void dump ( std::ostream& ss ) const { ss << mTriedge ; } ;

private :

  Triedge          mTriedge ;
  Trisegment_2_ptr mTrisegment ;
  Point_2          mP ;
  FT               mTime ;
} ;

template<class SSkel_, class Traits_>
class Edge_event_2 : public Event_2<SSkel_,Traits_>
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;

  typedef Event_2<SSkel,Traits> Base ;

  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;

  typedef typename Base::Type             Type ;
  typedef typename Base::Triedge          Triedge ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

public:

  Edge_event_2 ( Triedge const&          aTriedge
               , Trisegment_2_ptr const& aTrisegment 
               , Vertex_handle           aLSeed
               , Vertex_handle           aRSeed
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
  
  typedef typename Base::Type             Type ;
  typedef typename Base::Triedge          Triedge ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

public:

  Split_event_2 ( Triedge const&          aTriedge
                , Trisegment_2_ptr const& aTrisegment 
                , Vertex_handle           aSeed
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

  typedef typename Base::Type             Type ;
  typedef typename Base::Triedge          Triedge ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

public:

  Pseudo_split_event_2 ( Triedge const&          aTriedge
                       , Trisegment_2_ptr const& aTrisegment 
                       , Vertex_handle           aSeed0
                       , Vertex_handle           aSeed1
                       , bool                    aOppositeIs0
                       )
    :
      Base(aTriedge,aTrisegment)
    , mSeed0(aSeed0)
    , mSeed1(aSeed1)
    , mOppositeIs0(aOppositeIs0)
  {}

  virtual Type type() const { return this->cPseudoSplitEvent ; }

  virtual Vertex_handle seed0() const { return mSeed0 ; }
  virtual Vertex_handle seed1() const { return mSeed1 ; }

  bool opposite_node_is_seed_0() const { return mOppositeIs0 ; }
  
  bool is_at_source_vertex() const { return opposite_node_is_seed_0() ; }
  
  Vertex_handle opposite_seed() const { return opposite_node_is_seed_0() ? seed0() : seed1() ; }
  
private :

  virtual void dump ( std::ostream& ss ) const
  {
    this->Base::dump(ss);
    
    ss << " ("
       << "Seed0=" << mSeed0->id() << (  mOppositeIs0 ? " {Opp} " : " " ) 
       << "Seed1=" << mSeed1->id() << ( !mOppositeIs0 ? " {Opp}"  : "" ) 
       << ")" ;
  }

private :

  Vertex_handle mSeed0  ;
  Vertex_handle mSeed1 ;
  bool          mOppositeIs0 ;
} ;

}

} // end namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_EVENTS_2_H //
// EOF //
