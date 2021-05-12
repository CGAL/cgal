// Copyright (c) 2011 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Laurent Rineau


#ifndef CGAL_INTERNAL_STATIC_FILTERS_IS_DEGENERATE_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_IS_DEGENERATE_3_H

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {


template < typename K_base, typename SFK >
class Is_degenerate_3
  : public K_base::Is_degenerate_3
{
  typedef typename K_base::Ray_3     Ray_3;
  typedef typename K_base::Segment_3 Segment_3;
  typedef typename K_base::Plane_3 Plane_3;
  typedef typename K_base::Is_degenerate_3 Base;
  typedef typename K_base::Construct_source_3 Construct_source_3;
  typedef typename K_base::Construct_target_3 Construct_target_3;
  typedef typename K_base::Construct_second_point_3 Construct_second_point_3;
  typedef typename SFK::Equal_3 Equal_3;

public:

  typedef typename Base::result_type  result_type;

  using Base::operator();

  result_type
  operator()(const Segment_3& s) const
  {
    return Equal_3()(Construct_source_3()(s), Construct_target_3()(s));
  }


  result_type
  operator()(const Ray_3& r) const
  {
    return Equal_3()(Construct_source_3()(r), Construct_second_point_3()(r));
  }

  result_type
  operator()(const Plane_3& p) const
  {
    CGAL_BRANCH_PROFILER(std::string("semi-static attempts/calls to   : ") +
                         std::string(CGAL_PRETTY_FUNCTION), tmp);

    Get_approx<Plane_3> get_approx; // Identity functor for all planes
                                    // but lazy planes
    double a, b, c;

    if (fit_in_double(get_approx(p).a(), a) && fit_in_double(get_approx(p).b(), b) &&
        fit_in_double(get_approx(p).c(), c) )
    {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      return a == 0 && b == 0 && c == 0;
    }
    return Base::operator()(p);
  }


}; // end class Is_degenerate_3

} // end namespace Static_filters_predicates

} // end namespace internal

} // end namespace CGAL

#endif  // CGAL_INTERNAL_STATIC_FILTERS_IS_DEGENERATE_3_H
