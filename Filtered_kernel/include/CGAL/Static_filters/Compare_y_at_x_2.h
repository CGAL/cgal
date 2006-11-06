// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
//
// Author(s)     : Andreas Meyer

#ifndef CGAL_STATIC_FILTERS_COMPARE_Y_AT_X_2_H 
#define CGAL_STATIC_FILTERS_COMPARE_Y_AT_X_2_H

#include <CGAL/Profile_counter.h>
#include <CGAL/Static_filter_error.h>

CGAL_BEGIN_NAMESPACE

template < typename K_base, typename Kernel >
class SF_Compare_y_at_x_2
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
    
    CGAL_kernel_precondition(p.x() >= min(s.source().x(), s.target().x()) &&
                             p.x() <= max(s.source().x(), s.target().x()));
    
    if( Kernel().less_x_2_object()( s.source(), s.target() ) )
      return enum_cast<Comparison_result>(Kernel().orientation_2_object()
                                          (p, s.source(), s.target() ));
    else if ( Kernel().less_x_2_object()( s.target(), s.source() ) )
      return enum_cast<Comparison_result>(Kernel().orientation_2_object()
                                          (p, s.target(), s.source() ));
    else {
      if( p.y() < min(s.target().y(), s.source().y()) )
        return SMALLER;
      if( max(s.target().y(), s.source().y()) < p.y() )
        return LARGER;
      return EQUAL;
    }
  }
};

CGAL_END_NAMESPACE

#endif
