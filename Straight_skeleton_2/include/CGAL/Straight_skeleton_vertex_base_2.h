// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_VERTEX_BASE_2_H
#define CGAL_STRAIGHT_SKELETON_VERTEX_BASE_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_halfedge_base_2.h>
#include <CGAL/Trisegment_2.h>

#include <CGAL/circulator.h>
#include <CGAL/Origin.h>
#include <CGAL/use.h>

#include <boost/iterator/iterator_facade.hpp>

#include <limits>

namespace CGAL {

template < class Refs, class P, class N >
class Straight_skeleton_vertex_base_base_2
{
  enum Flags { IsSplitBit = 0x01, HasInfiniteTimeBit = 0x02 } ;

protected :

  class Halfedge_circulator_around_vertex_access_policy
  {
  public:
    template<class Impl>
    static typename Impl::reference access ( Impl* aImpl )
    {
      return aImpl->mHandle;
    }
  } ;

  class Halfedge_circulator_across_incident_faces_access_policy
  {
  public:
    template<class Impl>
    static typename Impl::reference access ( Impl* aImpl )
    {
      return aImpl->mHandle->face()->halfedge();
    }
  } ;

  template<class HalfedgeHandle, class AccessPolicy >
  class Halfedge_circulator_base
    : public boost::iterator_facade< Halfedge_circulator_base< HalfedgeHandle, AccessPolicy >
                                   ,HalfedgeHandle
                                   ,Bidirectional_circulator_tag
                                   ,HalfedgeHandle
                                  >
  {
    public:

      typedef HalfedgeHandle               value_type ;
      typedef HalfedgeHandle               reference ;
      typedef std::size_t                  size_type ;
      typedef Bidirectional_circulator_tag iterator_category ;

      Halfedge_circulator_base () : mHandle() {}

      explicit Halfedge_circulator_base ( value_type aHandle ) : mHandle(aHandle) {}

      template < class OtherHalfedgeHandle, class OtherAccessPolicy >
      Halfedge_circulator_base
        ( Halfedge_circulator_base<OtherHalfedgeHandle,OtherAccessPolicy> const& aOther )
        : mHandle(aOther.mHandle) {}

      bool operator==( std::nullptr_t p ) const
      {
        CGAL_USE(p);
        CGAL_assertion( p == nullptr );
        HalfedgeHandle null ;
        return mHandle == null ;
      }

      bool operator!=( std::nullptr_t p ) const { return !(*this == p); }

    private :

      typedef Halfedge_circulator_base<HalfedgeHandle,AccessPolicy> Self ;

      friend class boost::iterator_core_access ;

      template < class OtherHalfedgeHandle, class OtherAccessPolicy >
      bool equal( Halfedge_circulator_base<OtherHalfedgeHandle,OtherAccessPolicy> const& aOther ) const
      {
        return mHandle == aOther.mHandle;
      }

      void increment() { mHandle = mHandle->opposite()->prev(); }

      void decrement() { mHandle = mHandle->next()->opposite() ; }

      reference dereference() const { return AccessPolicy::access(const_cast<Self*>(this)) ; }

    private :

      friend class Halfedge_circulator_around_vertex_access_policy ;
      friend class Halfedge_circulator_across_incident_faces_access_policy ;

      value_type mHandle ;
  } ;

public:

  typedef Straight_skeleton_vertex_base_base_2<Refs, P, N>  Base ;

  typedef P Point_2;
  typedef N FT ;

  typedef Refs HalfedgeDS;

  typedef Tag_true Supports_vertex_halfedge;
  typedef Tag_true Supports_vertex_point;

  typedef typename Refs::Vertex_handle         Vertex_handle;
  typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
  typedef typename Refs::Halfedge_handle       Halfedge_handle;
  typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Refs::Face_handle           Face_handle;
  typedef typename Refs::Face_const_handle     Face_const_handle;
  typedef typename Refs::Halfedge              Halfedge;
  typedef typename Refs::Face                  Face;

  typedef Halfedge_circulator_base< Halfedge_const_handle
                                   ,Halfedge_circulator_around_vertex_access_policy
                                  >
            Halfedge_around_vertex_const_circulator ;

  typedef Halfedge_circulator_base< Halfedge_handle
                                   ,Halfedge_circulator_around_vertex_access_policy
                                  >
            Halfedge_around_vertex_circulator ;

  typedef Halfedge_circulator_base< Halfedge_const_handle
                                   ,Halfedge_circulator_across_incident_faces_access_policy
                                  >
            Defining_contour_halfedges_const_circulator ;

  typedef Halfedge_circulator_base< Halfedge_handle
                                   ,Halfedge_circulator_across_incident_faces_access_policy
                                  >
            Defining_contour_halfedges_circulator ;

  typedef CGAL_SS_i::Triedge<Halfedge_handle> Triedge ;

  typedef typename CGAL::Kernel_traits<P>::type K ;
  typedef CGAL_SS_i::Segment_2_with_ID<K> Segment_2 ;
  typedef CGAL_SS_i::Segment_2_with_ID<K> Segment_2_with_ID ; // for BOOST_MPL_HAS_XXX_TRAIT_DEF
  typedef CGAL::Trisegment_2<K, Segment_2_with_ID> Trisegment_2 ;
  typedef CGAL::Trisegment_2_ptr<Trisegment_2> Trisegment_2_ptr;

public:

  Straight_skeleton_vertex_base_base_2() : mID(-1), mTime(0.0), mFlags(0) {}

  // Infinite vertex
  Straight_skeleton_vertex_base_base_2 ( int aID )
    :
      mID   (aID)
    , mP    (ORIGIN)
    , mTime ((std::numeric_limits<double>::max)())
    , mFlags(HasInfiniteTimeBit)
  {
  }

  // Contour vertex
  Straight_skeleton_vertex_base_base_2 ( int aID, Point_2 const& aP )
    :
      mID   (aID)
    , mP    (aP)
    , mTime (0.0)
    , mFlags(0)
  {
  }

  // Skeleton vertex, corresponding to a split or edge event.
  Straight_skeleton_vertex_base_base_2 ( int aID, Point_2 const& aP, FT aTime, bool aIsSplit, bool aHasInfiniteTime )
    :
      mID   ( aID )
    , mP    ( aP )
    , mTime ( aTime )
    , mFlags( ( aIsSplit ? IsSplitBit : 0 ) | ( aHasInfiniteTime ? HasInfiniteTimeBit : 0 ) )
 {
 }

public:

  int id() const { return mID ; }

  FT time() const { return mTime ; }

  bool has_infinite_time() const { return ( mFlags & HasInfiniteTimeBit ) == HasInfiniteTimeBit ; }

  bool is_split() const { return ( mFlags & IsSplitBit ) == IsSplitBit ; }

  Halfedge_const_handle primary_bisector() const { return halfedge()->next(); }

  Halfedge_handle primary_bisector() { return halfedge()->next(); }

  Halfedge_around_vertex_const_circulator halfedge_around_vertex_begin() const
  {
    return Halfedge_around_vertex_const_circulator(halfedge());
  }

  Halfedge_around_vertex_circulator halfedge_around_vertex_begin()
  {
    return Halfedge_around_vertex_circulator(halfedge());
  }

  Defining_contour_halfedges_const_circulator defining_contour_halfedges_begin() const
  {
    return Defining_contour_halfedges_const_circulator(halfedge());
  }

  Defining_contour_halfedges_circulator defining_contour_halfedges_begin()
  {
    return Defining_contour_halfedges_circulator(halfedge());
  }

  bool is_skeleton() const { return  halfedge()->is_bisector() ; }
  bool is_contour () const { return !halfedge()->is_bisector() ; }

  const Point_2& point() const { return mP; }

  Halfedge_handle       halfedge()       { return mHE; }
  Halfedge_const_handle halfedge() const { return mHE; }

  void set_halfedge( Halfedge_handle aHE)  { mHE = aHE; }

  // Store a pointer to the trisegment, which also includes its potential children.
  // This is done to keep in memory the history of each node as to be able to
  // recompute its geometric position and time during offset polygon construction.
  //
  // Note: the trisegment stored was constructed in the straight skeleton builder.
  // When FinishUp() is called, multinodes are processed but as nodes are merged,
  // the trisegments of these nodes are *not* updated. Thus, the combinatorial trees
  // of these trisegments will become incoherent with the straight skeleton, but
  // that's OK because it is still valid to compute purely geometrical information
  // such as the node position and its time, which is all that is required for offset tracing.
  Trisegment_2_ptr const& trisegment() const { return mTrisegment ; }
  void set_trisegment( Trisegment_2_ptr const& aTrisegment ) { mTrisegment = aTrisegment ; }

public :

  void reset_id__internal__    ( int aID ) { mID = aID ; }
  void reset_point__internal__ ( Point_2 const& aP ) { mP = aP ; }

private:

  int              mID ;
  Halfedge_handle  mHE ;
  Point_2          mP;
  FT               mTime ;
  unsigned char    mFlags ;
  Trisegment_2_ptr mTrisegment ;
};

template < class Refs, class P, class N >
class Straight_skeleton_vertex_base_2 : public Straight_skeleton_vertex_base_base_2<Refs,P,N>
{
public:

  typedef P Point_2;
  typedef N FT ;

  typedef typename Refs::Vertex_handle         Vertex_handle;
  typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
  typedef typename Refs::Halfedge_handle       Halfedge_handle;
  typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Refs::Face_handle           Face_handle;
  typedef typename Refs::Face_const_handle     Face_const_handle;
  typedef typename Refs::Halfedge              Halfedge;
  typedef typename Refs::Face                  Face;

  typedef Straight_skeleton_vertex_base_base_2<Refs,P,N> Base ;

  typedef typename Base::Triedge               Triedge ;
  typedef typename Base::Trisegment_2_ptr      Trisegment_2_ptr ;

  Straight_skeleton_vertex_base_2() {}

  Straight_skeleton_vertex_base_2 ( int aID ) : Base(aID) {}

  Straight_skeleton_vertex_base_2 ( int aID, Point_2 const& aP ) : Base(aID,aP) {}

  Straight_skeleton_vertex_base_2 ( int aID, Point_2 const& aP, FT aTime, bool aIsSplit, bool aHasInfiniteTime )
    : Base(aID, aP, aTime, aIsSplit, aHasInfiniteTime)
  {}
} ;

} // namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_VERTEX_BASE_2_H //
