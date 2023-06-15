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
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_CACHES_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_CACHES_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>

#include <CGAL/assertions.h>

#include <vector>

namespace CGAL {
namespace CGAL_SS_i {

template <class Info>
struct FPU_checker;

template<class NT>
class Rational;

template <class Info>
struct No_cache
{
  bool IsCached ( std::size_t ) const { return false; }

  Info Get ( std::size_t ) const
  {
    CGAL_error();
    return Info();
  }

  void Set ( std::size_t, Info const& ) { }
  void Reset ( std::size_t ) { }
};

//TODO: call reserve, but how? #input vertices + n*m estimation?
template <class Info>
struct Info_cache
{
  std::vector<Info> mValues ;
  std::vector<bool> mAlreadyComputed ;

  bool IsCached ( std::size_t i ) const
  {
    return ( (mAlreadyComputed.size() > i) && mAlreadyComputed[i] ) ;
  }

  Info const& Get(std::size_t i) const
  {
    CGAL_precondition ( IsCached(i) ) ;
    CGAL_precondition ( FPU_checker<Info>::is_valid() ) ;
    return mValues[i] ;
  }

  void Set ( std::size_t i, Info const& aValue)
  {
    CGAL_precondition ( FPU_checker<Info>::is_valid() ) ;
    if (mValues.size() <= i )
    {
      mValues.resize(i+1) ;
      mAlreadyComputed.resize(i+1, false) ;
    }

    mAlreadyComputed[i] = true ;
    mValues[i] = aValue ;
  }

  void Reset ( std::size_t i )
  {
    if ( IsCached(i) ) // needed if approx info is set but not exact info
      mAlreadyComputed[i] = false ;
  }
};

template <typename K>
using Coeff_cache = Info_cache<std::optional<typename K::Line_2> > ;

template <typename K>
using Time_cache = Info_cache<std::optional<Rational<typename K::FT> > > ;

template <typename K>
using Point_cache = Info_cache<std::optional<typename K::Point_2> > ;

template <typename K>
struct Caches
{
  void Reset ( std::size_t i )
  {
    // !WARNING! subtlety here:
    // - The coefficient caches are attached to *contour edge IDs*.
    //   There is never any reason to reset coefficient caches because the geometry does not change.
    // - The time and point caches are attached to *trisegment IDs*
    //   There is a reason to reset trisegment caches to not waste memory when a trisegment is ignored.
    //
    // mCoeff_cache.Reset(i) ;

    mTime_cache.Reset(i) ;
    mPoint_cache.Reset(i) ;
  }

  Coeff_cache<K> mCoeff_cache;
  Time_cache<K> mTime_cache;
  Point_cache<K> mPoint_cache;
};

template <typename K>
struct No_caches
{
  void Reset ( std::size_t ) { }

  No_cache<std::optional<typename K::Line_2> > mCoeff_cache;
  No_cache<std::optional<Rational<typename K::FT> > > mTime_cache;
  No_cache<std::optional<typename K::Point_2> > mPoint_cache;
};

} // namespace CGAL_SS_i
} // namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_CACHES_H
