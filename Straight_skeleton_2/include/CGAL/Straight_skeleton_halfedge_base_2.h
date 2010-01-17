// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
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
#ifndef CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H
#define CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H 1


CGAL_BEGIN_NAMESPACE

template < class Refs, class FT_>
class Straight_skeleton_halfedge_base_base_2
{
protected:

  enum Flags { IsBisectorBit = 0x02, SlopeBitmask = 0x0C, PositiveSlope = 0x04, NegativeSlope = 0x08, ZeroSlope = 0x0C } ;

public:

  typedef Refs     HalfedgeDS;
  typedef Tag_true Supports_halfedge_prev;
  typedef Tag_true Supports_halfedge_vertex;
  typedef Tag_true Supports_halfedge_face;
  
  typedef typename Refs::Vertex_handle         Vertex_handle;
  typedef typename Refs::Vertex_const_handle   Vertex_const_handle;
  typedef typename Refs::Halfedge_handle       Halfedge_handle;
  typedef typename Refs::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Refs::Face_handle           Face_handle;
  typedef typename Refs::Face_const_handle     Face_const_handle;
  typedef typename Refs::Vertex                Vertex;
  typedef typename Refs::Face                  Face;
  
  typedef FT_ FT ;
  
  typedef Straight_skeleton_halfedge_base_base_2<Refs,FT> Base_base ;
  
  
protected:

  Straight_skeleton_halfedge_base_base_2( int aID, unsigned char aFlags, FT aW ) :  mF(Face_handle()), mID(aID), mFlags(aFlags), mW(aW) {}
  
public:

  int id() const { return mID ; }
  
  bool is_bisector() const { return test_bit(IsBisectorBit) ; }
  bool is_boundary() const { return ! is_bisector(); }

  bool is_inner_bisector() const
  {
    return is_bisector() && !vertex()->is_contour() && !opposite()->vertex()->is_contour();
  }

  bool has_null_segment() const { return !handle_assigned(vertex()) || vertex()->has_null_point() ; }
  
  bool has_infinite_time() const { return !handle_assigned(vertex()) || vertex()->has_infinite_time() ; }
  
  Halfedge_const_handle defining_contour_edge() const { return handle_assigned(face()) ? face()->halfedge() : Halfedge_const_handle() ; }
  Halfedge_handle       defining_contour_edge()       { return handle_assigned(face()) ? face()->halfedge() : Halfedge_handle() ; }

  Halfedge_handle       opposite()       { return mOpp;}
  Halfedge_const_handle opposite() const { return mOpp;}
  Halfedge_handle       next    ()       { return mNxt;}
  Halfedge_const_handle next    () const { return mNxt;}
  Halfedge_handle       prev    ()       { return mPrv; }
  Halfedge_const_handle prev    () const { return mPrv; }
  Vertex_handle         vertex  ()       { return mV; }
  Vertex_const_handle   vertex  () const { return mV; }
  Face_handle           face    ()       { return mF; }
  Face_const_handle     face    () const { return mF; }
  FT const&             weight  () const { return mW; }
  
  bool has_positive_slope() const { return test_bits(SlopeBitmask, PositiveSlope) ; }
  bool has_negative_slope() const { return test_bits(SlopeBitmask, NegativeSlope) ; }
  bool has_zero_slope    () const { return test_bits(SlopeBitmask, ZeroSlope    ) ; }
  
  Sign slope() const { return has_positive_slope() ? POSITIVE : ( has_negative_slope() ? NEGATIVE : ZERO ) ; }

  void set_opposite( Halfedge_handle h) { mOpp = h; }
  void set_next    ( Halfedge_handle h) { mNxt = h; }
  void set_prev    ( Halfedge_handle h) { mPrv = h; }
  void set_vertex  ( Vertex_handle   w) { mV   = w; }
  void set_face    ( Face_handle     g) { mF   = g; }
  void set_weight  ( FT const&       w) { mW   = w ; }
 
  void set_slope( Sign aSlope ) { set_bits(SlopeBitmask, slope_bits(aSlope) ) ; }

  void reset_id ( int aID ) { mID = aID ; }

  bool has_no_incident_face() const { return !handle_assigned(face()); }

private:
  
  void set_is_bisector( bool aOn ) { set_bits( IsBisectorBit, aOn ? IsBisectorBit : 0 ) ; }
  
  void set_bits  ( unsigned char aBitmask, unsigned char aBits ) { mFlags = ( mFlags & ~aBitmask ) | ( aBitmask & aBits ) ;  }
  
  bool test_bits ( unsigned char aBitmask, unsigned char aBits ) const { return ( mFlags & aBitmask ) == aBits ; }
  
  bool test_bit  ( unsigned char aBit ) const { return test_bits(aBit,aBit); }
  
  static unsigned char slope_bits( Sign aSlope ) { return aSlope == POSITIVE ? PositiveSlope : ( aSlope == NEGATIVE ? NegativeSlope : ZeroSlope ) ; }
  
private:

  Halfedge_handle mOpp;
  Halfedge_handle mNxt;
  Halfedge_handle mPrv;
  Vertex_handle   mV;
  Face_handle     mF;
  FT              mW;
  int             mID ;
  unsigned char   mFlags ;
};

template < class Refs, class FT >
class Straight_skeleton_halfedge_base_2 : public Straight_skeleton_halfedge_base_base_2<Refs,FT>
{
public:

  typedef typename Refs::Vertex_handle   Vertex_handle;
  typedef typename Refs::Halfedge_handle Halfedge_handle;
  typedef typename Refs::Face_handle     Face_handle;
  
  typedef Straight_skeleton_halfedge_base_base_2<Refs,FT> Base_base ;
  typedef Straight_skeleton_halfedge_base_2<Refs,FT>      Base ;     

  Straight_skeleton_halfedge_base_2( int aID = -1, unsigned char aFlags = 0, FT aW = 1.0 ) : Base_base(aID, aFlags, aW) {}
  
  static Straight_skeleton_halfedge_base_2 NewBoundary( int aID )
  {
    return Straight_skeleton_halfedge_base_2(aID, 0);
  }
  
  static Straight_skeleton_halfedge_base_2 NewBisector( int aID )
  {
    return Straight_skeleton_halfedge_base_2(aID, Base::IsBisectorBit );
  }
  
private:

  void set_opposite( Halfedge_handle h )  { Base_base::opposite(h)   ; }
  void set_next    ( Halfedge_handle h )  { Base_base::set_next(h)   ; }
  void set_prev    ( Halfedge_handle h )  { Base_base::set_prev(h)   ; }
  void set_vertex  ( Vertex_handle   v )  { Base_base::set_vertex(v) ; }
  void set_face    ( Face_handle     f )  { Base_base::set_face(f)   ; }
  void set_slope   ( Sign            s )  { Base_base::set_slope(s)  ; }
  void set_weight  ( FT const&       w )  { Base_base::set_weight(w) ; }
  void reset_id    ( int             i )  { Base_base::reset_id(i)   ; }

} ;

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_HALFEDGE_BASE_2_H //
// EOF //

