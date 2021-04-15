// Copyright (c) 2017  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_STATICALLY_FILTERED_TRAITS_2_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_STATICALLY_FILTERED_TRAITS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/internal/Static_filters/Periodic_2_side_of_oriented_circle_2.h>
#include <CGAL/internal/Periodic_2_Delaunay_triangulation_filtered_traits_2.h>

namespace CGAL {

template<class K_,
         class Off_ = typename CGAL::Periodic_2_offset_2>
class Periodic_2_Delaunay_triangulation_statically_filtered_traits_2
  : public Periodic_2_Delaunay_triangulation_filtered_traits_base_2<K_, Off_>
{
  typedef Periodic_2_Delaunay_triangulation_statically_filtered_traits_2<K_, Off_> Self;
  typedef Periodic_2_Delaunay_triangulation_filtered_traits_base_2<K_, Off_>       Base;

public:
  typedef K_                                                                       Kernel;
  typedef typename Kernel::Iso_rectangle_2                                         Iso_rectangle_2;

  Periodic_2_Delaunay_triangulation_statically_filtered_traits_2(const Iso_rectangle_2& domain,
                                                                 const Kernel& k)
    : Base(domain, k)
  { }

  typedef internal::Static_filters_predicates::Periodic_2_side_of_oriented_circle_2<
            Self, typename Base::Side_of_oriented_circle_2> Side_of_oriented_circle_2;

  Side_of_oriented_circle_2  side_of_oriented_circle_2_object() const {
    return Side_of_oriented_circle_2(&this->_domain,
                                     this->Base::side_of_oriented_circle_2_object());
  }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_STATICALLY_FILTERED_TRAITS_2_H
