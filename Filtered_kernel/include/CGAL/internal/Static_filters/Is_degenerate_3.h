// Copyright (c) 2011 GeometryFactory (France)
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
  typedef typename K_base::Is_degenerate_3 Base;
  typedef typename K_base::Construct_source_3 Construct_source_3;
  typedef typename K_base::Construct_target_3 Construct_target_3;
  typedef typename K_base::Construct_second_point_3 Construct_second_point_3;
  typedef typename SFK::Equal_3 Equal_3;

public:

  typedef typename Base::result_type  result_type;


#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else // CGAL_CFG_MATCHING_BUG_6
  template <typename T>
  result_type
  operator()(const T& t) const
  {
    return Base()(t);
  }
#endif // end CGAL_CFG_MATCHING_BUG_6


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

}; // end class Is_degenerate_3

} // end namespace Static_filters_predicates

} // end namespace internal

} // end namespace CGAL

#endif  // CGAL_INTERNAL_STATIC_FILTERS_IS_DEGENERATE_3_H
