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

template<class SSkel_, class Traits_>
class Straight_skeleton_builder_event_2 : public Ref_counted_base
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;
  
public:
 
  typedef Straight_skeleton_builder_event_2<SSkel,Traits> Self ;

  typedef boost::intrusive_ptr<Self> SelfPtr ;

  typedef typename Traits::Point_2          Point_2 ;
  typedef typename Traits::FT               FT ;
  typedef typename Traits::Sorted_triedge_2 Sorted_triedge_2 ;

  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;

  enum Type { cEdgeEvent, cSplitEvent, cPseudoSplitEvent } ;

public:

  Straight_skeleton_builder_event_2 (  Halfedge_handle         aBorderA
                                     , Halfedge_handle         aBorderB
                                     , Halfedge_handle         aBorderC
                                     , Sorted_triedge_2 const& aSorted 
                                    )
    :
     mBorderA(aBorderA)
    ,mBorderB(aBorderB)
    ,mBorderC(aBorderC)
    ,mSorted (aSorted)
    ,mExcluded(false)
  {}

  virtual ~ Straight_skeleton_builder_event_2() {}

  virtual Type type() const = 0 ;

  virtual Vertex_handle seed0() const = 0 ;
  virtual Vertex_handle seed1() const = 0 ;

  Halfedge_handle         border_a      () const { return mBorderA ; }
  Halfedge_handle         border_b      () const { return mBorderB ; }
  Halfedge_handle         border_c      () const { return mBorderC ; }
  Sorted_triedge_2 const& sorted_triedge() const { return mSorted  ; }
  Point_2 const&          point         () const { return mP       ; }
  FT                      time          () const { return mTime    ; }

  bool is_excluded() const { return mExcluded ; }
  void Exclude    ()       { mExcluded = true ; }

  void SetTimeAndPoint( FT aTime, Point_2 const& aP ) { mTime = aTime ; mP = aP ; }

  friend std::ostream& operator<< ( std::ostream& ss, Self const& e )
  {
    ss << "[" ;
    e.dump(ss);
    ss << " p=(" << e.point().x() << "," << e.point().y() << ") t=" << e.time() << "] " 
       << triedge_collinearity_to_string(e.sorted_triedge().collinearity()) ;
    return ss ;
  }

protected :

  virtual void dump ( std::ostream& ss ) const
  {
    ss << "{E" << mBorderA->id() << ",E" << mBorderB->id() << ",E" << mBorderC->id() << "}" ;
  } ;

private :

  Halfedge_handle  mBorderA ;
  Halfedge_handle  mBorderB ;
  Halfedge_handle  mBorderC ;
  Sorted_triedge_2 mSorted ;
  Point_2          mP ;
  FT               mTime ;
  bool             mExcluded ;
} ;

template<class SSkel_, class Traits_>
class Straight_skeleton_builder_edge_event_2 : public Straight_skeleton_builder_event_2<SSkel_,Traits_>
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;

  typedef Straight_skeleton_builder_event_2<SSkel,Traits> Base ;

  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;

  typedef typename Base::Type Type ;
  
  typedef typename Traits::Sorted_triedge_2 Sorted_triedge_2 ;

public:

  Straight_skeleton_builder_edge_event_2 (  Halfedge_handle         aBorderA
                                          , Halfedge_handle         aBorderB
                                          , Halfedge_handle         aBorderC
                                          , Sorted_triedge_2 const& aSorted 
                                          , Vertex_handle           aLSeed
                                          , Vertex_handle           aRSeed
                                          )
    :
      Base(aBorderA,aBorderB,aBorderC,aSorted)
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
class Straight_skeleton_builder_split_event_2 : public Straight_skeleton_builder_event_2<SSkel_,Traits_>
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;

  typedef Straight_skeleton_builder_event_2<SSkel,Traits> Base ;

  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;
  
  typedef typename Base::Type Type ;

  typedef typename Traits::Sorted_triedge_2 Sorted_triedge_2 ;

public:

  Straight_skeleton_builder_split_event_2 (  Halfedge_handle         aBorderA
                                           , Halfedge_handle         aBorderB
                                           , Halfedge_handle         aBorderC
                                           , Sorted_triedge_2 const& aSorted 
                                           , Vertex_handle           aSeed
                                         )
    :
      Base(aBorderA,aBorderB,aBorderC,aSorted)
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
    ss << " (Seed=" << mSeed->id() << " OppBorder=" << this->border_c()->id() << ')' ;
  }

private :

  Vertex_handle mSeed ;
  Vertex_handle mOppR ;
} ;

template<class SSkel_, class Traits_>
class Straight_skeleton_builder_pseudo_split_event_2 : public Straight_skeleton_builder_event_2<SSkel_,Traits_>
{
  typedef SSkel_  SSkel ;
  typedef Traits_ Traits ;
  
  typedef Straight_skeleton_builder_event_2<SSkel,Traits> Base ;

  typedef typename SSkel::Halfedge_handle Halfedge_handle ;
  typedef typename SSkel::Vertex_handle   Vertex_handle ;

  typedef typename Base::Type Type ;
  
  typedef typename Traits::Sorted_triedge_2 Sorted_triedge_2 ;

public:

  Straight_skeleton_builder_pseudo_split_event_2 ( Halfedge_handle         aBorderA
                                                 , Halfedge_handle         aBorderB
                                                 , Halfedge_handle         aBorderC
                                                 , Sorted_triedge_2 const& aSorted 
                                                 , Vertex_handle           aSeed
                                                 , Vertex_handle           aOppositeNode
                                                 )
    :
      Base(aBorderA,aBorderB,aBorderC,aSorted)
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


CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_EVENTS_2_H //
// EOF //
