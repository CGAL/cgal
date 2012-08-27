// Copyright (c) 2004   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_REGULAR_TRIANGULATION_FILTERED_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_FILTERED_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Filtered_predicate.h>

namespace CGAL { 

// The Weighted_converter is parametrized by a usual kernel converter,
// and adds the conversions for the Weighted_point.
template < typename Converter >
struct Weighted_converter_2
  : Converter
{
  typedef typename Converter::Source_kernel Source_kernel;
  typedef typename Converter::Target_kernel Target_kernel;

  typedef typename Regular_triangulation_euclidean_traits_base_2<Source_kernel>
                   ::Weighted_point_2  Source_wp;

  typedef typename Regular_triangulation_euclidean_traits_base_2<Target_kernel>
                   ::Weighted_point_2  Target_wp;

  Target_wp
  operator()(const Source_wp &wp) const
  {
    return Target_wp(Converter::operator()(wp.point()),
                     Converter::operator()(wp.weight()));
  }
};

// The argument is supposed to be a Filtered_kernel like kernel.
template < typename K >
class Regular_triangulation_filtered_traits_2
  : public Regular_triangulation_euclidean_traits_base_2<K>
{
  // Exact traits is based on the exact kernel.
  typedef Regular_triangulation_euclidean_traits_2<typename K::Exact_kernel>
                                                   Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Regular_triangulation_euclidean_traits_2<typename K::Approximate_kernel>
                                                   Filtering_traits;

  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

public:

  typedef K               Kernel;

  typedef Filtered_predicate<
            typename Exact_traits::Power_test_2,
            typename Filtering_traits::Power_test_2,
            Weighted_converter_2<C2E>,
            Weighted_converter_2<C2F> >  Power_test_2;

  typedef Filtered_predicate<
    typename Exact_traits::Compare_power_distance_2,
    typename Filtering_traits::Compare_power_distance_2,
    Weighted_converter_2<C2E>,
    Weighted_converter_2<C2F> >  Compare_power_distance_2;

  Power_test_2 power_test_2_object() const
  { return Power_test_2();}

  Compare_power_distance_2 compare_power_distance_2_object() const
  { return Compare_power_distance_2(); }

  // The following are inherited since they are constructions :
  // Construct_weighted_circumcenter_2
  // Construct_radical_axis_2
};

} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_FILTERED_TRAITS_2_H
