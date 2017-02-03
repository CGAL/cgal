// Copyright (c) 1999-2004,2006-2009,2014-2015   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Aymeric Pellé <Aymeric.Pelle@sophia.inria.fr>
#ifndef CGAL_REGULAR_TRAITS_WITH_OFFSETS_ADAPTOR_H
#define CGAL_REGULAR_TRAITS_WITH_OFFSETS_ADAPTOR_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/Traits_with_offsets_adaptor.h>

namespace CGAL {

template < class PRTT, class Functor_ >
class Regular_traits_with_offsets_adaptor
  : public Traits_with_offsets_adaptor<PRTT, Functor_>
{
  typedef PRTT                                              PRTraits;
  typedef Functor_                                          Functor;
  typedef Traits_with_offsets_adaptor<PRTraits, Functor_>   Base;

  typedef typename PRTraits::FT                             FT;
  typedef typename PRTraits::Point_3                        Point_3;
  typedef typename PRTraits::Weighted_point_3               Weighted_point_3;
  typedef typename PRTraits::Offset                         Offset;
  typedef typename PRTraits::Construct_weighted_point_3     Construct_weighted_point_3;

public:
  typedef typename PRTraits::Iso_cuboid_3                   Iso_cuboid_3;
  typedef typename Functor::result_type                     result_type;

public:
  Regular_traits_with_offsets_adaptor (const Iso_cuboid_3 * dom)
    : Base(dom)
  { }

public:
  using Base::operator();
  using Base::pp;

  Weighted_point_3 pp(const Weighted_point_3 &p, const Offset &o) const {
    return Construct_weighted_point_3(*(this->_domain))(p, o);
  }

public:
  // with offset ---------------------------------------------------------------
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Offset& o0, const Offset& o1) const
  {
    return Functor()(pp(p0, o0), pp(p1, o1));
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2) const
  {
    return Functor()(pp(p0, o0), pp(p1, o1), pp(p2, o2));
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2, const Offset& o3) const
  {
    return Functor()(pp(p0, o0), pp(p1, o1), pp(p2, o2), pp(p3, o3));
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const Weighted_point_3& p4,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2, const Offset& o3,
                          const Offset& o4) const
  {
    return Functor()(pp(p0, o0), pp(p1, o1), pp(p2, o2), pp(p3, o3), pp(p4, o4));
  }

  // for `Compare_power_distance_3`
  result_type operator() (const Point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2) const
  {
    return Functor()(pp(p0, o0), pp(p1, o1), pp(p2, o2));
  }

  // for `Compare_weighted_squared_radius_3`
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2, const Offset& o3,
                          const FT w) const
  {
    return Functor()(pp(p0, o0), pp(p1, o1), pp(p2, o2), pp(p3, o3), w);
  }

  // without offset ------------------------------------------------------------
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1) const
  {
    return Functor()(p0, p1);
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2) const
  {
    return Functor()(p0, p1, p2);
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3) const
  {
    return Functor()(p0, p1, p2, p3);
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const Weighted_point_3& p4) const
  {
    return Functor()(p0, p1, p2, p3, p4);
  }

  // for `Compare_power_distance_3`
  result_type operator() (const Point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2) const
  {
    return Functor()(p0, p1, p2);
  }

  // for `Compare_weighted_squared_radius_3`
  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const FT w) const
  {
    return Functor()(p0, p1, p2, p3, w);
  }
};

}  // namespace CGAL

#endif /* CGAL_REGULAR_TRAITS_WITH_OFFSETS_ADAPTOR_H */
