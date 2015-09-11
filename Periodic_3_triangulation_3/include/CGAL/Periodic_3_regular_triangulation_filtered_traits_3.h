// Copyright (c) 1999-2004,2006-2009,2013-2015   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@inria.fr>
//                 Aymeric Pelle <Aymeric.Pelle@sophia.inria.fr>


#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_FILTERED_TRAITS_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_FILTERED_TRAITS_3_H

#include <string>
#include <CGAL/basic.h>
#include <CGAL/config.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/Periodic_3_triangulation_filtered_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/internal/Regular_triangulation_filtered_traits_3.h>

namespace CGAL {

template < typename K, typename Off >
class Periodic_3_regular_triangulation_filtered_traits_base_3
  : public Periodic_3_regular_triangulation_traits_base_3<K, Off>
{
  typedef Periodic_3_regular_triangulation_traits_3<K, Off> Base;

  // Exact traits is based on the exact kernel.
  typedef Periodic_3_regular_triangulation_traits_3<typename K::Exact_traits,
                                            Off>
                                                   Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef Periodic_3_regular_triangulation_traits_3<typename K::Filtering_traits,
                                            Off>
                                                   Filtering_traits;
private:
  typedef typename K::Kernel::C2E C2E;
  typedef typename K::Kernel::C2F C2F;

  typedef typename C2E::Target_kernel::Iso_cuboid_3 Exact_iso_cuboid_3;
  typedef typename C2F::Target_kernel::Iso_cuboid_3 Approximate_iso_cuboid_3;
 
public:
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

  void set_domain(const Iso_cuboid_3& domain) {
    C2E c2e;
    C2F c2f;
    this->_domain = domain;
    this->_domain_e = c2e(this->_domain);
    this->_domain_f = c2f(this->_domain);
  }

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Compare_xyz_3,
            typename Filtering_traits::Compare_xyz_3,
            Offset_converter_3<Weighted_converter_3<C2E> >,
            Offset_converter_3<Weighted_converter_3<C2F> > >  Compare_xyz_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Coplanar_orientation_3,
            typename Filtering_traits::Coplanar_orientation_3,
            Offset_converter_3<Weighted_converter_3<C2E> >,
            Offset_converter_3<Weighted_converter_3<C2F> > >  Coplanar_orientation_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Orientation_3,
            typename Filtering_traits::Orientation_3,
            Offset_converter_3<Weighted_converter_3<C2E> >,
            Offset_converter_3<Weighted_converter_3<C2F> > >  Orientation_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Power_test_3,
            typename Filtering_traits::Power_test_3,
            Offset_converter_3<Weighted_converter_3<C2E> >,
            Offset_converter_3<Weighted_converter_3<C2F> > >  Power_test_3;

  typedef Filtered_periodic_predicate<
            typename Exact_traits::Compare_power_distance_3,
            typename Filtering_traits::Compare_power_distance_3,
            Offset_converter_3<Weighted_converter_3<C2E> >,
            Offset_converter_3<Weighted_converter_3<C2F> > >  Compare_power_distance_3;

  Compare_xyz_3 compare_xyz_3_object() const
  { return Compare_xyz_3(&_domain_e,&_domain_f);}

  Coplanar_orientation_3 coplanar_orientation_3_object() const
  { return Coplanar_orientation_3(&_domain_e,&_domain_f); }

  Orientation_3 orientation_3_object() const
  { return Orientation_3(&_domain_e,&_domain_f);}

  Power_test_3 power_test_3_object () const
  {
    return Power_test_3(&_domain_e,&_domain_f);
  }

  Compare_power_distance_3 compare_power_distance_3_object () const
  {
    return Compare_power_distance_3(&_domain_e,&_domain_f);
  }

 protected:
  Exact_iso_cuboid_3 _domain_e;
  Approximate_iso_cuboid_3 _domain_f;
};

} //namespace CGAL

//#include <CGAL/Periodic_3_regular_triangulation_statically_filtered_traits_3.h>

namespace CGAL {

template < typename K, typename Off = typename CGAL::Periodic_3_offset_3 >
class Periodic_3_regular_triangulation_filtered_traits_3
    : public Periodic_3_regular_triangulation_filtered_traits_base_3<K, Off>
{
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_REGULAR_TRIANGULATION_FILTERED_TRAITS_3_H
