// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Pm_straight_exact_traits.C
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ============================================================================
#ifndef CGAL_PM_STRAIGHT_EXACT_TRAITS_C
#define CGAL_PM_STRAIGHT_EXACT_TRAITS_C

#ifndef CGAL_PM_STRAIGHT_EXACT_TRAITS_H
#include <CGAL/Pm_straight_exact_traits.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class R_>
const typename Pm_straight_exact_traits<R_>::Bounding_box 
Pm_straight_exact_traits<R_>::unbounded_box()
{
  static const typename Pm_straight_exact_traits<R_>::Bounding_box 
    unbounded_box_;
  return unbounded_box_;
}

template <class R_>
const typename Pm_straight_exact_traits<R_>::Bounding_box 
Pm_straight_exact_traits<R_>::unbounded_box_;
// unbounded_box_ initialized to default Bounding_box

CGAL_END_NAMESPACE

#endif
