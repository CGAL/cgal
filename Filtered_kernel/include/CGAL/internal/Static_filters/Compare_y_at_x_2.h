// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Meyer

#ifndef CGAL_INTERNAL_STATIC_FILTERS_COMPARE_Y_AT_X_2_H
#define CGAL_INTERNAL_STATIC_FILTERS_COMPARE_Y_AT_X_2_H

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <CGAL/algorithm.h>

namespace CGAL { namespace internal { namespace Static_filters_predicates {

template < typename K_base, typename Kernel >
class Compare_y_at_x_2
  : public K_base::Compare_y_at_x_2
{
  typedef typename K_base::Point_2          Point_2;
  typedef typename K_base::Segment_2        Segment_2;
  typedef typename K_base::FT               FT;
  typedef typename K_base::Compare_y_at_x_2 Base;

public:

  using Base::operator();

  Comparison_result
  operator()( const Point_2& p, const Segment_2& s ) const {
    // compares the y-coordinates of p and the vertical projection of p on s.
    // Precondition : p is in the x-range of s.

    typename Kernel::Less_x_2 less_x = Kernel().less_x_2_object();
    typename Kernel::Less_y_2 less_y = Kernel().less_y_2_object();
    typename Kernel::Orientation_2 orientation = Kernel().orientation_2_object();

    CGAL_kernel_precondition( are_ordered(s.source(), p, s.target(), less_x) );

    if( less_x( s.source(), s.target() ) )
      return orientation(p, s.source(), s.target());
    else if ( less_x( s.target(), s.source() ) )
      return orientation(p, s.target(), s.source());
    else {
      if( less_y(p, s.source()) && less_y(p, s.target()) )
        return SMALLER;
      if( less_y(s.source(), p) && less_y(s.target(), p) )
        return LARGER;
      return EQUAL;
    }
  }
};

} } } // namespace CGAL::internal::Static_filters_predicates

#endif // CGAL_INTERNAL_STATIC_FILTERS_COMPARE_Y_AT_X_2_H
