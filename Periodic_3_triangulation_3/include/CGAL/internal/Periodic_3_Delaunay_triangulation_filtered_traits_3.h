// Copyright (c) 2004,2006-2009, 2017  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>


#ifndef CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_FILTERED_TRAITS_3_H
#define CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_FILTERED_TRAITS_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <string>
#include <CGAL/basic.h>
#include <CGAL/config.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/internal/Periodic_3_triangulation_filtered_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>

namespace CGAL {

// The first template item is supposed to be a Filtered_kernel-like kernel.
template < typename K, typename Off >
class Periodic_3_Delaunay_triangulation_filtered_traits_base_3
  : public Periodic_3_Delaunay_triangulation_traits_base_3<K, Off>
{
  typedef Periodic_3_Delaunay_triangulation_traits_base_3<K, Off> Base;

  // Exact traits is based on the exact kernel.
  typedef Periodic_3_Delaunay_triangulation_traits_3<typename K::Exact_kernel,
                                                     Off>      Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Periodic_3_Delaunay_triangulation_traits_3<typename K::Approximate_kernel,
                                                     Off>      Filtering_traits;
private:
  typedef typename K::C2E C2E;
  typedef typename K::C2F C2F;

public:
  typedef Filtered_periodic_predicate<
            typename Exact_traits::Side_of_oriented_sphere_3,
            typename Filtering_traits::Side_of_oriented_sphere_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Side_of_oriented_sphere_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Compare_distance_3,
            typename Filtering_traits::Compare_distance_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Compare_distance_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Coplanar_orientation_3,
            typename Filtering_traits::Coplanar_orientation_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Coplanar_orientation_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Coplanar_side_of_bounded_circle_3,
            typename Filtering_traits::Coplanar_side_of_bounded_circle_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Coplanar_side_of_bounded_circle_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Side_of_bounded_sphere_3,
            typename Filtering_traits::Side_of_bounded_sphere_3,
            Offset_converter_3<C2E>,
            Offset_converter_3<C2F> >  Side_of_bounded_sphere_3;

  Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object() const
  { return Side_of_oriented_sphere_3(&this->_domain_e,&this->_domain_f); }

  Compare_distance_3 compare_distance_3_object() const
  { return Compare_distance_3(&this->_domain_e,&this->_domain_f); }

  Coplanar_orientation_3 coplanar_orientation_3_object() const
  { return Coplanar_orientation_3(&this->_domain_e,&this->_domain_f); }

  Coplanar_side_of_bounded_circle_3 coplanar_side_of_bounded_circle_3_object() const
  { return Coplanar_side_of_bounded_circle_3(&this->_domain_e,&this->_domain_f); }

  Side_of_bounded_sphere_3 side_of_bounded_sphere_3_object() const
  { return Side_of_bounded_sphere_3(&this->_domain_e,&this->_domain_f); }
};

} //namespace CGAL

#include <CGAL/internal/Periodic_3_Delaunay_triangulation_statically_filtered_traits_3.h>

namespace CGAL {

template < typename K,
           typename Off = typename CGAL::Periodic_3_offset_3,
           bool Has_static_filters = internal::Has_static_filters<K>::value >
class Periodic_3_Delaunay_triangulation_filtered_traits_3
  : public Periodic_3_Delaunay_triangulation_statically_filtered_traits_3<
      Periodic_3_Delaunay_triangulation_filtered_traits_base_3<K, Off> >
{ };

template < typename K_, typename Off>
class Periodic_3_Delaunay_triangulation_filtered_traits_3<K_, Off, false>
  :  public Periodic_3_Delaunay_triangulation_filtered_traits_base_3<K_, Off>
{ };

} //namespace CGAL

#endif // CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_FILTERED_TRAITS_3_H
