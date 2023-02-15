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


#ifndef CGAL_INTERNAL_STATIC_FILTERS_COPLANAR_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_COPLANAR_3_H

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {


template < typename K_base, typename SFK >
class Coplanar_3
  : public K_base::Coplanar_3
{
  typedef typename K_base::Point_3 Point_3;
  typedef typename K_base::Coplanar_3 Base;
  typedef typename SFK::Orientation_3 Orientation_3;

public:

  typedef typename Base::result_type  result_type;



  result_type
  operator()(const Point_3& p,const Point_3& q, const Point_3& r, const Point_3& s) const
  {
    return Orientation_3()(p,q,r,s) == COPLANAR;
  }



}; // end class Coplanar_3

} // end namespace Static_filters_predicates

} // end namespace internal

} // end namespace CGAL

#endif  // CGAL_INTERNAL_STATIC_FILTERS_COPLANAR_3_H
