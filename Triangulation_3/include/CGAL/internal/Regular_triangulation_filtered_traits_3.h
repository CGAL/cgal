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
// Author(s)     : Sylvain Pion

#ifndef CGAL_INTERNAL_REGULAR_TRIANGULATION_FILTERED_TRAITS_3_H
#define CGAL_INTERNAL_REGULAR_TRIANGULATION_FILTERED_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/internal/Static_filters/Regular_triangulation_static_filters_traits_3.h>

namespace CGAL{

  
// The Weighted_converter is parametrized by a usual kernel converter,
// and adds the conversions for the Weighted_point.
template < typename Converter, 
           typename Source_traits= Regular_triangulation_euclidean_traits_base_3<typename Converter::Source_kernel>,
           typename Target_traits= Regular_triangulation_euclidean_traits_base_3<typename Converter::Target_kernel> 
         >
struct Weighted_converter_3
  : public Converter
{
  typedef typename Converter::Source_kernel Source_kernel;
  typedef typename Converter::Target_kernel Target_kernel;

  typedef typename Source_traits::Weighted_point_3  Source_wp;
  typedef typename Target_traits::Weighted_point_3  Target_wp;

  typedef typename Source_kernel::Point_3  Source_p;
  typedef typename Target_kernel::Point_3  Target_p;

#ifdef CGAL_CFG_MATCHING_BUG_7

  typedef typename Source_kernel::FT  Source_FT;
  typedef typename Target_kernel::FT  Target_FT;

  
  // needed for weighted Alpha shapes
 Target_FT
  operator()(const Source_FT &p) const
  {
    return Converter::operator()(p);
  }
  
#else 
  using Converter::operator();
#endif 

  // Needed for MSVC 2005/2008 to avoid a matching ambiguity  
  Target_p
  operator()(const Source_p &p) const
  {
    return Converter::operator()(p);
  }
  
  Target_wp
  operator()(const Source_wp &wp) const
  {
    return Target_wp(Converter::operator()(wp.point()),
                     Converter::operator()(wp.weight()));
  }

};


namespace internal{

// The argument is supposed to be a Filtered_kernel like kernel.
template < typename K>
class Regular_triangulation_filtered_traits_base_3
  : public Regular_triangulation_euclidean_traits_base_3<K>
{
  // Exact traits is based on the exact kernel.
  typedef Regular_triangulation_euclidean_traits_3<typename K::Exact_kernel>
                                                   Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Regular_triangulation_euclidean_traits_3<typename K::Approximate_kernel>
                                                   Filtering_traits;

  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

public:

  typedef K               Kernel;

  typedef Filtered_predicate<
            typename Exact_traits::Power_test_3,
            typename Filtering_traits::Power_test_3,
            Weighted_converter_3<C2E>,
            Weighted_converter_3<C2F> >  Power_test_3;

  typedef Filtered_predicate<
            typename Exact_traits::Compare_power_distance_3,
            typename Filtering_traits::Compare_power_distance_3,
            Weighted_converter_3<C2E>,
            Weighted_converter_3<C2F> >  Compare_power_distance_3;

  typedef Filtered_predicate<
            typename Exact_traits::In_smallest_orthogonal_sphere_3,
            typename Filtering_traits::In_smallest_orthogonal_sphere_3,
            Weighted_converter_3<C2E>,
            Weighted_converter_3<C2F> >  In_smallest_orthogonal_sphere_3;

  typedef Filtered_predicate<
            typename Exact_traits::Side_of_bounded_orthogonal_sphere_3,
            typename Filtering_traits::Side_of_bounded_orthogonal_sphere_3,
            Weighted_converter_3<C2E>,
            Weighted_converter_3<C2F> >  Side_of_bounded_orthogonal_sphere_3;

  typedef Filtered_predicate<
            typename Exact_traits::Does_simplex_intersect_dual_support_3,
            typename Filtering_traits::Does_simplex_intersect_dual_support_3,
            Weighted_converter_3<C2E>,
            Weighted_converter_3<C2F> >  Does_simplex_intersect_dual_support_3;

 typedef Filtered_predicate<
            typename Exact_traits::Compare_weighted_squared_radius_3,
            typename Filtering_traits::Compare_weighted_squared_radius_3,
            Weighted_converter_3<C2E>,
            Weighted_converter_3<C2F> >  Compare_weighted_squared_radius_3;

  enum { Has_filtered_predicates = true };
  
  Power_test_3 power_test_3_object() const
  { return Power_test_3();}

  Compare_power_distance_3 compare_power_distance_3_object() const
  { return Compare_power_distance_3();}

  In_smallest_orthogonal_sphere_3
  in_smallest_orthogonal_sphere_3_object() const
  { return In_smallest_orthogonal_sphere_3(); }

  Side_of_bounded_orthogonal_sphere_3
  side_of_bounded_orthogonal_sphere_3_object() const
  { return Side_of_bounded_orthogonal_sphere_3(); }

  Does_simplex_intersect_dual_support_3
  does_simplex_intersect_dual_support_3_object() const
  { return Does_simplex_intersect_dual_support_3(); }

  Compare_weighted_squared_radius_3
  compare_weighted_squared_radius_3_object() const
  { return Compare_weighted_squared_radius_3();  }  
  // The following are inherited since they are constructions :
  // Construct_weighted_circumcenter_3
  // Compute_squared_radius_smallest_orthogonal_sphere_3
  // Compute_power_product_3
};

template < typename K,bool UseStaticFilters = K::Has_static_filters>
class Regular_triangulation_filtered_traits_3
  : public Regular_triangulation_filtered_traits_base_3<K>
{
public:
  enum { Has_static_filters = false };
};

template < typename K>
class Regular_triangulation_filtered_traits_3<K,true>
    : public internal::Regular_triangulation_static_filters_traits_3< Regular_triangulation_filtered_traits_base_3<K> >
{
public:  
	enum { Has_static_filters = true };
};

	

} } //namespace CGAL::internal

#endif // CGAL_INTERNAL_REGULAR_TRIANGULATION_FILTERED_TRAITS_3_H
