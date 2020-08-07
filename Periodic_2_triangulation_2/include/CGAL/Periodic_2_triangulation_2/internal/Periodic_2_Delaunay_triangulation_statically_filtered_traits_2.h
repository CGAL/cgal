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

#include <CGAL/Periodic_2_triangulation_2/internal/Static_filters/Periodic_2_side_of_oriented_circle_2.h>
#include <CGAL/Periodic_2_triangulation_2/internal/Periodic_2_Delaunay_triangulation_filtered_traits_2.h>

#include <CGAL/Default.h>

namespace CGAL {

template<class Kernel_,
         class Offset_ = typename CGAL::Periodic_2_offset_2,
         class Domain_ = typename Kernel_::Iso_rectangle_2>
class Periodic_2_Delaunay_triangulation_statically_filtered_traits_2
  : public Periodic_2_Delaunay_triangulation_filtered_traits_base_2<Kernel_, Offset_, Domain_>
{
  typedef Periodic_2_Delaunay_triangulation_statically_filtered_traits_2<Kernel_, Offset_, Domain_> Self;
  typedef Periodic_2_Delaunay_triangulation_filtered_traits_base_2<Kernel_, Offset_, Domain_>       Base;

public:
  typedef Kernel_                                                                        Kernel;
  typedef Domain_                                                                        Domain;

  Periodic_2_Delaunay_triangulation_statically_filtered_traits_2(const Domain& domain,
                                                                 const Kernel& k = Kernel())
    : Base(domain, k)
  {
    std::cout << "create 'Periodic_2_Delaunay_triangulation_statically_filtered_traits_2'" << std::endl;
  }

  typedef internal::Static_filters_predicates::Periodic_2_side_of_oriented_circle_2<
            Self, typename Base::Side_of_oriented_circle_2> Side_of_oriented_circle_2;

  Side_of_oriented_circle_2  side_of_oriented_circle_2_object() const {
    return Side_of_oriented_circle_2(&this->domain, this->Base::side_of_oriented_circle_2_object());
  }
};

// @todo static filters for lattices
template<class Kernel_, class Offset_>
class Periodic_2_Delaunay_triangulation_statically_filtered_traits_2<Kernel_, Offset_, Lattice_2<Kernel_> >
  : public Periodic_2_Delaunay_triangulation_filtered_traits_base_2<Kernel_, Offset_, Lattice_2<Kernel_> >
{
public:
  typedef Kernel_                                                                        Kernel;
  typedef Lattice_2<Kernel_>                                                             Domain;

private:
  typedef Periodic_2_Delaunay_triangulation_filtered_traits_base_2<Kernel_, Offset_, Domain> Base;

public:
  Periodic_2_Delaunay_triangulation_statically_filtered_traits_2(const Domain& domain = P2T2::internal::get_default_domain<Kernel, Domain>(),
                                                                 const Kernel& k = Kernel())
    : Base(domain, k)
  { }
};

} // namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_STATICALLY_FILTERED_TRAITS_2_H
