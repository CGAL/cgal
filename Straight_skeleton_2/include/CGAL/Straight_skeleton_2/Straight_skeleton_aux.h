// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_AUX_H
#define CGAL_STRAIGHT_SKELETON_AUX_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Straight_skeleton_2/assertions.h>
#include <CGAL/Straight_skeleton_2/debug.h>
#include <CGAL/Straight_skeleton_2/test.h>

#include <CGAL/Polygon_with_holes_2.h>

#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/or.hpp>
#include <optional>

#include <iostream>
#include <type_traits>

namespace CGAL {

namespace CGAL_SS_i {

template<class K>
struct Has_inexact_constructions
{
  typedef typename K::FT FT ;

  typedef std::conditional_t< std::is_same_v<FT,double> || std::is_same_v<FT,Interval_nt_advanced>
                              , Tag_true
                              , Tag_false
                              > type ;
} ;

template <class K>
struct Segment_2_with_ID
  : public K::Segment_2
{
  typedef typename K::Segment_2 Base;
  typedef typename K::Point_2 Point_2;

  std::size_t id() const { return mID; }

public:
  Segment_2_with_ID() : Base(), mID(-1) { }
  Segment_2_with_ID(Base const& aS) : Base(aS), mID(-1) { }
  Segment_2_with_ID(Base const& aS, const std::size_t aID) : Base(aS), mID(aID) { }
  Segment_2_with_ID(Point_2 const& aP, Point_2 const& aQ, const std::size_t aID) : Base(aP, aQ), mID(aID) { }

public:
  std::size_t mID;
};

//
// This record encapsulates the defining contour halfedges for a node (both contour and skeleton)
//
template<class Handle_>
class Triedge
{
public:

  typedef Handle_ Handle ;

  typedef Triedge<Handle> Self ;

  Triedge() {}

  // Contour nodes (input polygon vertices) have only 2 defining contour edges
  Triedge ( Handle aE0, Handle aE1 )
  {
    mE[0] = aE0 ;
    mE[1] = aE1 ;
    // mE[2] gets default constructed, i.e., "null".
  }

  // Skeleton nodes (offset polygon vertices) have 3 defining contour edges
  Triedge ( Handle aE0, Handle aE1 , Handle aE2 )
  {
    mE[0] = aE0 ;
    mE[1] = aE1 ;
    mE[2] = aE2 ;
  }

  Handle e( unsigned idx ) const { CGAL_assertion(idx<3); return mE[idx]; }

  Handle e0() const { return e(0); }
  Handle e1() const { return e(1); }
  Handle e2() const { return e(2); }

  bool is_valid() const { return handle_assigned(e0()) && handle_assigned(e1()); }

  bool is_contour () const { return !handle_assigned(e2()) ; }
  bool is_skeleton() const { return  handle_assigned(e2()) ; }

  bool is_contour_terminal() const { return e0() == e1() ; }

  bool is_skeleton_terminal() const { return e0() == e1() || e1() == e2() ; }

  // Returns true if the triedges store the same 3 halfedges (in any order)
  friend bool operator == ( Self const& x, Self const& y )
  {
    return x.number_of_unique_edges() == y.number_of_unique_edges() && CountInCommon(x,y) == x.number_of_unique_edges() ;
  }

  friend bool operator != ( Self const& x, Self const& y ) { return !(x==y) ; }

  friend Self operator & ( Self const& x, Self const& y )
  {
    return Self(x.e0(), x.e1(), ( x.e0() == y.e0() || x.e1() == y.e0() ) ? y.e1() : y.e0()  ) ;
  }

  friend std::ostream& operator<< ( std::ostream& ss, Self const& t )
  {
    ss << "{E" ;
    insert_handle_id(ss,t.e0())  ;
    ss << ",E" ;
    insert_handle_id(ss,t.e1()) ;
    ss << ",E" ;
    insert_handle_id(ss,t.e2()) ;
    ss << "}" ;
    return ss ;
  }

private:

  static void insert_handle_id( std::ostream& ss, Handle aH )
  {
    if ( handle_assigned(aH) )
         ss << aH->id() ;
    else ss << "#" ;
  }

  // returns 1 if aE is one of the halfedges stored in this triedge, 0 otherwise.
  int contains ( Handle aE ) const
  {
    return aE == e0() || aE == e1() || aE == e2() ? 1 : 0 ;
  }

  int number_of_unique_edges() const
  {
    return is_contour() ? ( is_contour_terminal() ? 1 : 2 ) : ( is_skeleton_terminal() ? 2 : 3 ) ;
  }

  // Returns the number of common halfedges in the two triedges x and y
  static int CountInCommon( Self const& x, Self const& y )
  {
    Handle lE[3];

    int lC = 1 ;

    lE[0] = y.e0();

    if ( y.e0() != y.e1() )
      lE[lC++] = y.e1();

    if ( y.e0() != y.e2() && y.e1() != y.e2() )
       lE[lC++] = y.e2();

    return x.contains(lE[0]) + x.contains(lE[1]) + ( lC > 2 ? x.contains(lE[2]) : 0 ) ;
  }

  Handle mE[3];
} ;

template<class U, class V>
struct Is_same_type { typedef Tag_false type ; } ;

template<class U>
struct Is_same_type<U,U> { typedef Tag_true type ; } ;

// to distinguish between SequenceContainers of points, and GeneralPolygonWithHoles_2
BOOST_MPL_HAS_XXX_TRAIT_DEF(Hole_const_iterator)

// The return type of create_(weighted)_interior/exterior_skeleton_and_offset_polygons_2:
// - if polygon input is a model of 'GeneralPolygonWithHoles_2', the return type
//   should be the internal (hole-less) polygon type GeneralPolygonWithHoles_2::Polygon_2
// - if polygon input is just a sequence container of points (e.g. Polygon_2), then the same type
//   is expected in output
template <typename Polygon, typename OfK,
          bool has_holes = CGAL_SS_i::has_Hole_const_iterator<Polygon>::value>
struct Default_return_polygon_type // Polygon type supports holes
{
  typedef typename std::conditional<std::is_same<
                                      typename Kernel_traits<typename boost::range_value<
                                        typename Polygon::Polygon_2>::type>::Kernel,
                                      OfK>::value,
                                    typename Polygon::Polygon_2, // correct kernel
                                    CGAL::Polygon_2<OfK> /*incorrect kernel*/ >::type type;
};

template <typename Polygon, typename OfK>
struct Default_return_polygon_type<Polygon, OfK, false> // Polygon type does NOT support holes
{
  typedef typename std::conditional<std::is_same<
                                      typename Kernel_traits<typename boost::range_value<Polygon>::type>::Kernel,
                                      OfK>::value,
                                    Polygon, // correct kernel
                                    CGAL::Polygon_2<OfK> /*incorrect kernel*/ >::type type;
};

// The return type of create_interior/exterior_skeleton_and_offset_polygons_with_holes_2:
// - if polygon input is a model of 'GeneralPolygonWithHoles_2', the return type should be the same
// - if polygon input is just a sequence container of points (e.g. Polygon_2), then use
//   General_polygon_with_holes_2<Polygon>
template <typename Polygon, typename OfK,
          bool has_holes = CGAL_SS_i::has_Hole_const_iterator<Polygon>::value>
struct Default_return_polygon_with_holes_type // Polygon type supports holes
{
  typedef typename std::conditional<std::is_same<
                                      typename Kernel_traits<typename boost::range_value<
                                        typename Polygon::Polygon_2>::type>::Kernel,
                                      OfK>::value,
                                    Polygon, // correct kernel
                                    CGAL::Polygon_with_holes_2<OfK> /*incorrect kernel*/ >::type type;
};

template <typename Polygon, typename OfK>
struct Default_return_polygon_with_holes_type<Polygon, OfK, false> // Polygon type does NOT support holes
{
  // Maybe on paper the `conditional<true>` should be `General_polygon_with_holes_2<Polygon>`...
  typedef typename std::conditional<std::is_same<
                                      typename Kernel_traits<typename boost::range_value<Polygon>::type>::Kernel,
                                      OfK>::value,
                                    CGAL::Polygon_with_holes_2<OfK>, // correct kernel but no holes
                                    CGAL::Polygon_with_holes_2<OfK> /*incorrect kernel*/ >::type type;
};

} // namespace CGAL_SS_i

} // namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_AUX_H
