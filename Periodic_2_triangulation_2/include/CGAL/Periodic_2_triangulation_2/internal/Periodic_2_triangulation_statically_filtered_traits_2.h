// Copyright (c) 1997-2013, 2017 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>,
//                 Mael Rouxel-Labbé

// This class gathers optimized predicates written by hand, using
// a few steps of filtering.  It should work if the initial traits has
// cartesian coordinates which fit exactly in doubles.
//
// Purely static filters code has been removed, since it requires additional
// logic and is not plug'n play (requires users providing bounds).
// If it should be provided again, it should probably be separate.

#ifndef CGAL_PERIODIC_2_TRIANGULATION_STATICALLY_FILTERED_TRAITS_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_STATICALLY_FILTERED_TRAITS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_2/internal/Periodic_2_triangulation_filtered_traits_2.h>
#include <CGAL/Periodic_2_triangulation_2/internal/Static_filters/Periodic_2_orientation_2.h>

namespace CGAL {

template<class Kernel_,
         class Offset_ = typename CGAL::Periodic_2_offset_2,
         class Domain_ = typename Kernel_::Iso_rectangle_2>
class Periodic_2_triangulation_statically_filtered_traits_2
  : public Periodic_2_triangulation_filtered_traits_base_2<Kernel_, Offset_, Domain_>
{
  typedef Periodic_2_triangulation_statically_filtered_traits_2<Kernel_, Offset_, Domain_>  Self;
  typedef Periodic_2_triangulation_filtered_traits_base_2<Kernel_, Offset_, Domain_>        Base;

public:
  typedef Kernel_                                                                           Kernel;
  typedef Domain_                                                                           Domain;

  Periodic_2_triangulation_statically_filtered_traits_2(const Domain& domain = P2T2::internal::get_default_domain<Kernel, Domain>(),
                                                        const Kernel& k = Kernel())
    : Base(domain, k)
  {
    std::cout << "create 'Periodic_2_triangulation_statically_filtered_traits_2'" << std::endl;
  }

  typedef internal::Static_filters_predicates::Periodic_2_orientation_2<
            Self, typename Base::Orientation_2> Orientation_2;

  Orientation_2 orientation_2_object() const
  {
    return Orientation_2(&this->domain, this->Base::orientation_2_object());
  }
};

template<class K, class O>
class Periodic_2_triangulation_statically_filtered_traits_2<K, O, Lattice_2<K> >
  : public Periodic_2_triangulation_filtered_traits_base_2<K, O, Lattice_2<K> >
{
public:
  typedef K                                                                                 Kernel;
  typedef Lattice_2<K>                                                                      Domain;

private:
  typedef Periodic_2_triangulation_filtered_traits_base_2<K, O, Domain>                     Base;

public:
  Periodic_2_triangulation_statically_filtered_traits_2(const Domain& domain = P2T2::internal::get_default_domain<K, Domain>(),
                                                        const K& k = Kernel())
    : Base(domain, k)
  { }
};

} // namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_STATICALLY_FILTERED_TRAITS_2_H
