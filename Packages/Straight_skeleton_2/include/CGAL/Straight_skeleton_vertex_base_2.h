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
// file          : include/CGAL/Straight_skeleton_vertex_base_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================

#ifndef CGAL_STRAIGHT_SKELETON_VERTEX_BASE_2_H
#define CGAL_STRAIGHT_SKELETON_VERTEX_BASE_2_H 1

#ifndef CGAL_HALFEDGEDS_VERTEX_BASE_H
#include <CGAL/HalfedgeDS_vertex_base.h>
#endif

#ifndef CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H
#include <CGAL/Straight_skeleton_halfedge_base_2.h>
#endif

#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif

#include <boost/iterator/iterator_facade.hpp>

CGAL_BEGIN_NAMESPACE

template < class Refs, class P, class N >
class Straight_skeleton_vertex_base_base_2 
  : public HalfedgeDS_vertex_base<Refs, Tag_true, P >
{

protected :
  
  template<class HalfedgeHandle, class AccessPolicy >
  class Halfedge_circulator_base
    : boost::iterator_facade< Halfedge_circulator_base< HalfedgeHandle, AccessPolicy >
                             ,HalfedgeHandle
                             ,Bidirectional_circulator_tag
                            > 
  {
    public:
      
      typedef HalfedgeHandle value_type ;
      
      Halfedge_circulator_base () : mH() {}
      
      explicit Halfedge_circulator_base ( value_type aHandle ) : mHandle(aHandle) {}
      
      template < class OtherHalfedgeHandle, class OtherAccessPolicy >
      Halfedge_circulator_base 
        ( Halfedge_circulator_base<OtherHalfedgeHandle,OtherAccessPolicy> const& aOther )
        : mHandle(aOther.mHandle) {}
      
    private :
    
      friend class boost::iterator_core_access ;

      template <class,class> friend class Halfedge_circulator_base;
     
      template < class OtherHalfedgeHandle, class OtherAccessPolicy >
      bool equal( Halfedge_circulator_base<OtherHalfedgeHandle,OtherAccessPolicy> const& aOther ) const 
      {
        return mHandle == aOther.mHandle;
      }

      void increment() { mHandle = mHandle->opposite()->next(); }      
      
      void decrement() { mHandle = mHandle->prev()->opposite() ; }
      
      value_type& dereference() const { return *AccessPolicy::access(mHandle) ; }
      
    private :
    
      value_type mHandle ;
  } ;
 
  class Halfedge_circulator_around_vertex_access_policy
  {
    template<class HalfedgeHandle>
    static HalfedgeHandle access ( HalfedgeHandle aHandle )
    {
      return aHandle;
    }
  } ;
  
  class Halfedge_circulator_across_incident_faces_access_policy
  {
    template<class HalfedgeHandle>
    static HalfedgeHandle access ( HalfedgeHandle aHandle )
    {
      return aHandle->face()->halfedge();
    }
  } ;
  
public:
 
  typedef HalfedgeDS_vertex_base<Refs, Tag_true, P> Base ;
  
  typedef P Point_2;
  typedef N FT ;

  typedef typename Refs::Vertex_handle         Vertex_handle ;
  typedef typename Refs::Vertex_const_handle   Vertex_const_handle ;
  typedef typename Refs::Halfedge_handle       Halfedge_handle ;
  typedef typename Refs::Halfedge_const_handle Halfedge_const_handle ;

  typedef Halfedge_circulator_base< Halfedge_const_handle const
                                   ,Halfedge_circulator_around_vertex_access_policy
                                  >
            Halfedge_around_vertex_const_circulator ;
  
  typedef Halfedge_circulator_base< Halfedge_handle
                                   ,Halfedge_circulator_around_vertex_access_policy
                                  >
            Halfedge_around_vertex_circulator ;
            
  typedef Halfedge_circulator_base< Halfedge_const_handle const
                                   ,Halfedge_circulator_across_incident_faces_access_policy                                  
                                  >
            Halfedge_across_incident_faces_const_circulator ;
  
  typedef Halfedge_circulator_base< Halfedge_handle
                                   ,Halfedge_circulator_across_incident_faces_access_policy
                                  >
            Halfedge_across_incident_faces_circulator ;
            
protected:
  
  Straight_skeleton_vertex_base_base_2() : mID(-1) {}

  Straight_skeleton_vertex_base_base_2 ( int aID, Point const& aP )
    :
      Base(aP)
    , mID(aID)
    , mTime(0.0)
  {} 
  
  Straight_skeleton_vertex_base_base_2 ( int aID, Point const& aP, FT aTime )
    :
      Base(aP)
    , mID(aID)
    , mTime(aTime)
 {}

  
public:

  int id() const { return mID ; }
  
  FT time() const { return mTime ; }
    
  Halfedge_around_vertex_const_circulator incident_edges_begin() const
  {
    return Halfedge_around_vertex_const_circulator(halfedge()); 
  }
    
  Halfedge_around_vertex_circulator incident_edges_begin()
  {
    return Halfedge_around_vertex_circulator(halfedge()); 
  }

  Halfedge_across_incident_faces_const_circulator defining_borders_begin() const
  {
    return Halfedge_across_incident_faces_const_circulator(halfedge());
  }    
  
  Halfedge_across_incident_faces_circulator defining_borders_begin()
  {
    return Halfedge_across_incident_faces_circulator(halfedge());
  }    
    
  bool is_inner () const { return  halfedge()->is_bisector() ; }
  bool is_border() const { return !halfedge()->is_bisector() ; }
  
private:

  int  mID ;
  FT   mTime ;
};

template < class Refs, class P, class N >
class Straight_skeleton_vertex_base_2 
  : public Straight_skeleton_vertex_base_base_2<Refs,P,N> 
{
  
public:
 
  typedef Straight_skeleton_vertex_base_base_2<Refs,P,N> Base ;
  
  typedef typename Base::Base            Base_base ;
  typedef typename Base::Point_2         Point_2;
  typedef typename Base::FT              FT;
  typedef typename Base::Halfedge_handle Halfedge_handle ;
  
public:
  
  Straight_skeleton_vertex_base_2() {}

  Straight_skeleton_vertex_base_2 ( int aID, Point const& aP )
    :
    Base(aID,aP)
  {} 
  
  Straight_skeleton_vertex_base_2 ( int aID, Point const& aP, FT aTime )
    :
    Base(aID,aP,aTime)    
 {}

protected:

  void set_halfedge( Halfedge_handle h ) { Base_base::set_halfedge(h) ; } 
};

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_VERTEX_BASE_2_H //
// EOF //
 
