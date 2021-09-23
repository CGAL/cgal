// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//

#ifndef CGAL_SLS_TRISEGMENT_H
#define CGAL_SLS_TRISEGMENT_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>

#include <CGAL/assertions.h>

#include <boost/intrusive_ptr.hpp>

#include <limits>
#include <iostream>

namespace CGAL {

enum Trisegment_collinearity
{
    TRISEGMENT_COLLINEARITY_NONE
  , TRISEGMENT_COLLINEARITY_01
  , TRISEGMENT_COLLINEARITY_12
  , TRISEGMENT_COLLINEARITY_02
  , TRISEGMENT_COLLINEARITY_ALL
} ;


inline char const* trisegment_collinearity_to_string( Trisegment_collinearity c )
{
  switch ( c )
  {
    case TRISEGMENT_COLLINEARITY_NONE : return "<>" ;
    case TRISEGMENT_COLLINEARITY_01   : return "<0,1>" ;
    case TRISEGMENT_COLLINEARITY_12   : return "<1,2>" ;
    case TRISEGMENT_COLLINEARITY_02   : return "<0,2>" ;
    case TRISEGMENT_COLLINEARITY_ALL  : return "<0,1,2>" ;
  }

  return "!!UNKNOWN COLLINEARITY!!" ;
}

namespace internal
{

template <>
struct Minmax_traits< Trisegment_collinearity >
{
  static const Trisegment_collinearity min = TRISEGMENT_COLLINEARITY_NONE;
  static const Trisegment_collinearity max = TRISEGMENT_COLLINEARITY_ALL;
};

}

template<class K, typename Segment>
class Trisegment_2
  : public Ref_counted_base
{
  typedef Trisegment_2<K, Segment> Self;

public:
  typedef boost::intrusive_ptr<Trisegment_2> Self_ptr ;

  Trisegment_2 ( Segment const&        aE0
               , Segment const&        aE1
               , Segment const&        aE2
               , Trisegment_collinearity aCollinearity
               , std::size_t aID
               )
    : mID(aID)
  {
    mCollinearity = aCollinearity ;

    mE[0] = aE0 ;
    mE[1] = aE1 ;
    mE[2] = aE2 ;

    switch ( mCollinearity )
    {
      case TRISEGMENT_COLLINEARITY_01:
        mCSIdx=0; mNCSIdx=2; break ;

      case TRISEGMENT_COLLINEARITY_12:
        mCSIdx=1; mNCSIdx=0; break ;

      case TRISEGMENT_COLLINEARITY_02:
        mCSIdx=0; mNCSIdx=1; break ;

      case TRISEGMENT_COLLINEARITY_ALL:
        mCSIdx = mNCSIdx = (std::numeric_limits<unsigned>::max)(); break ;

      case TRISEGMENT_COLLINEARITY_NONE:
        mCSIdx = mNCSIdx = (std::numeric_limits<unsigned>::max)(); break ;
    }
  }

  std::size_t& id() { return mID; }
  const std::size_t& id() const { return mID; }

  static Trisegment_2 null() { return Self_ptr() ; }

  Trisegment_collinearity collinearity() const { return mCollinearity ; }

  Segment const& e( unsigned idx ) const { CGAL_precondition(idx<3) ; return mE[idx] ; }

  Segment const& e0() const { return e(0) ; }
  Segment const& e1() const { return e(1) ; }
  Segment const& e2() const { return e(2) ; }

  // If 2 out of the 3 edges are collinear they can be reclassified as 1 collinear edge (any of the 2) and 1 non-collinear.
  // These methods returns the edges according to that classification.
  // PRECONDITION: Exactly 2 out of 3 edges are collinear
  Segment const& collinear_edge    () const { return e(mCSIdx) ; }
  Segment const& non_collinear_edge() const { return e(mNCSIdx) ; }
  Segment const& other_collinear_edge() const
  {
    switch ( mCollinearity )
    {
      case TRISEGMENT_COLLINEARITY_01:
        return e(1);
      case TRISEGMENT_COLLINEARITY_12:
        return e(2);
      case TRISEGMENT_COLLINEARITY_02:
        return e(2);
    }
  }

  Self_ptr child_l() const { return mChildL ; }
  Self_ptr child_r() const { return mChildR ; }
  Self_ptr child_t() const { return mChildT ; }

  void set_child_l( Self_ptr const& aChild ) { mChildL = aChild ; }
  void set_child_r( Self_ptr const& aChild ) { mChildR = aChild ; }
  void set_child_t( Self_ptr const& aChild ) { mChildT = aChild ; }

  enum SEED_ID { LEFT, RIGHT, THIRD } ;

  // Indicates which of the seeds is collinear for a normal collinearity case.
  // PRECONDITION: The collinearity is normal.
  SEED_ID degenerate_seed_id() const
  {
    Trisegment_collinearity c = collinearity();

    return c == TRISEGMENT_COLLINEARITY_01 ? LEFT : c == TRISEGMENT_COLLINEARITY_12 ? RIGHT : THIRD  ;
  }

  friend std::ostream& operator << ( std::ostream& os, Self const& aTrisegment )
  {
    return os << "[" << s2str(aTrisegment.e0())
              << " " << s2str(aTrisegment.e1())
              << " " << s2str(aTrisegment.e2())
              << " " << trisegment_collinearity_to_string(aTrisegment.collinearity())
              << "]";
  }

  friend std::ostream& operator << ( std::ostream& os, Self_ptr const& aPtr )
  {
    recursive_print(os,aPtr,0);
    return os ;
  }

  static void recursive_print ( std::ostream& os, Self_ptr const& aTriPtr, int aDepth )
  {
    os << "\n" ;

    for ( int i = 0 ; i < aDepth ; ++ i )
      os << "  " ;

    if ( aTriPtr )
    {
      os << *aTriPtr ;

      if ( aTriPtr->child_l() )
        recursive_print(os,aTriPtr->child_l(),aDepth+1);

      if ( aTriPtr->child_r() )
        recursive_print(os,aTriPtr->child_r(),aDepth+1);
    }
    else
    {
      os << "{null}" ;
    }
  }

private :
  std::size_t             mID;
  Segment                 mE[3];
  Trisegment_collinearity mCollinearity ;
  unsigned                mCSIdx, mNCSIdx ;

  Self_ptr mChildL ;
  Self_ptr mChildR ;

  // this is the potential child of e2-e0, if it exists. It is used only in the configuration
  // of e0 and e2 collinear as the common child gives where the bisector starts (as it is not
  // necessarily the middle of the gap between e2 and e0).
  Self_ptr mChildT ;
} ;

} // end namespace CGAL

#endif // CGAL_SLS_TRISEGMENT_H
