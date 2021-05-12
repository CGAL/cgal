// Copyright (c) 2017  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_3_POWER_TEST_3_H
#define CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_3_POWER_TEST_3_H

#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Static_filters/Static_filter_error.h>
#include <CGAL/internal/Static_filters/tools.h>

#include <CGAL/Periodic_3_offset_3.h>

#include <cmath>

namespace CGAL {

namespace internal {

namespace Static_filters_predicates {

template <class K, class Power_side_of_oriented_power_sphere_3_base>
class Periodic_3_power_side_of_oriented_power_sphere_3:
    public Power_side_of_oriented_power_sphere_3_base
{
  typedef Power_side_of_oriented_power_sphere_3_base           Base;

public:
  typedef K                                                    Kernel;

  // @todo
};

} // namespace Static_filters_predicates

} // namespace internal

} //namespace CGAL

#endif // CGAL_INTERNAL_STATIC_FILTERS_PERIODIC_3_POWER_TEST_3_H
