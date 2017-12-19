// Copyright (c) 2008 ETH Zurich (Switzerland)
// Copyright (c) 2008-2009 INRIA Sophia-Antipolis (France)
// Copyright (c) 2017 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Andreas Fabri, Laurent Rineau


#ifndef CGAL_INTERNAL_STATIC_FILTERS_DO_INTERSECT_2_H
#define CGAL_INTERNAL_STATIC_FILTERS_DO_INTERSECT_2_H

#include <iostream>

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {

template < typename K_base, typename SFK >
class Do_intersect_2
  : public K_base::Do_intersect_2
{
  typedef typename K_base::Point_2   Point_2;
  typedef typename K_base::Segment_2 Segment_2;
  typedef typename K_base::Do_intersect_2 Base;

  typedef K_base TA1;
  typedef SFK TA2;
public:

  typedef typename Base::result_type  result_type;


#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else // CGAL_CFG_MATCHING_BUG_6
  template <typename T1, typename T2>
  result_type
  operator()(const T1& t1, const T2& t2) const
  {
    return Base()(t1,t2);
  }
#endif // CGAL_CFG_MATCHING_BUG_6

  // The internal::do_intersect(..) function 
  // only performs orientation tests on the vertices
  // of the segment
  // By calling the do_intersect function with
  // the  statically filtered kernel we avoid
  // that doubles are put into Interval_nt
  // to get taken out again with fit_in_double
  result_type 
  operator()(const Segment_2 &s, const Segment_2& t) const
  {
    return internal::do_intersect(s,t, SFK());
  }

  result_type 
  operator()(const Point_2 &p, const Segment_2& t) const
  {
    return internal::do_intersect(p,t, SFK());
  }
  
  result_type 
  operator()(const Segment_2& t, const Point_2 &p) const
  {
    return internal::do_intersect(p,t, SFK());
  }
  
};
} // Static_filters_predicates
} // internal
} // CGAL
#endif
