 // Copyright (c) 2009  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_INTERNAL_TRIANGULATION_2_FILTERED_PROJECTION_TRAITS_3_H
#define CGAL_INTERNAL_TRIANGULATION_2_FILTERED_PROJECTION_TRAITS_3_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/Triangulation_2/internal/Triangulation_2_projection_traits_base_3.h>
#include <CGAL/Filtered_predicate_with_state.h>

namespace CGAL {

template < class Filtered_kernel >
class Triangulation_2_filtered_projection_traits_3
  : public Triangulation_2_projection_traits_base_3<Filtered_kernel>
{
  typedef Filtered_kernel K;
  typedef Triangulation_2_filtered_projection_traits_3<K> Self;
  typedef Triangulation_2_projection_traits_base_3<K> Base;

  typedef typename K::Exact_kernel Exact_kernel;
  typedef typename K::Approximate_kernel Approximate_kernel;
  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

public:
  typedef Triangulation_2_projection_traits_base_3<Exact_kernel> Exact_traits;
  typedef Triangulation_2_projection_traits_base_3<Approximate_kernel> Filtering_traits;

public:
  explicit Triangulation_2_filtered_projection_traits_3(const typename K::Vector_3& n)
    : Base(n)
  {
  }

#define CGAL_TRIANGULATION_2_PROJ_TRAITS_FILTER_PRED(P, Pf, ACCESSOR)    \
  typedef  Filtered_predicate_with_state< \
    typename Exact_traits::P, \
    typename Filtering_traits::P, \
    C2E, \
    C2F, \
    typename K::Vector_3> P;      \
  P Pf() const { \
    return P(this->ACCESSOR()); \
  }
  CGAL_TRIANGULATION_2_PROJ_TRAITS_FILTER_PRED(Orientation_2,
                                               orientation_2_object,
                                               normal)
  CGAL_TRIANGULATION_2_PROJ_TRAITS_FILTER_PRED(Side_of_oriented_circle_2,
                                               side_of_oriented_circle_2_object,
                                               normal)
  CGAL_TRIANGULATION_2_PROJ_TRAITS_FILTER_PRED(Less_x_2,
                                               less_x_2_object,
                                               base1)
  CGAL_TRIANGULATION_2_PROJ_TRAITS_FILTER_PRED(Less_y_2,
                                               less_y_2_object,
                                               base2)
  CGAL_TRIANGULATION_2_PROJ_TRAITS_FILTER_PRED(Compare_x_2,
                                               compare_x_2_object,
                                               base1)
  CGAL_TRIANGULATION_2_PROJ_TRAITS_FILTER_PRED(Compare_y_2,
                                               compare_y_2_object,
                                               base2)
}; // end class Triangulation_2_projection_traits_base_3<Filtered_kernel>

} // end namespace CGAL


#endif // CGAL_INTERNAL_TRIANGULATION_2_FILTERED_PROJECTION_TRAITS_3_H
