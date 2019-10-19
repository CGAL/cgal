// Copyright (c) 1999-2004,2006-2009,2014-2015,2017  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Aymeric Pellé <Aymeric.Pelle@sophia.inria.fr>
//                 Mael Rouxel-Labbé
#ifndef CGAL_FUNCTOR_WITH_OFFSET_WEIGHTED_POINTS_ADAPTOR_3_H
#define CGAL_FUNCTOR_WITH_OFFSET_WEIGHTED_POINTS_ADAPTOR_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/internal/Functor_with_offset_points_adaptor_3.h>

namespace CGAL {

template < class K_, class Functor_ >
class Functor_with_offset_weighted_points_adaptor_3
  : public Functor_with_offset_points_adaptor_3<K_, Functor_>
{
  typedef Functor_with_offset_points_adaptor_3<K_, Functor_>  Base;

  typedef K_                                                Kernel;
  typedef Functor_                                          Functor;

  typedef typename Kernel::FT                               FT;
  typedef typename Kernel::Point_3                          Point_3;
  typedef typename Kernel::Weighted_point_3                 Weighted_point_3;
  typedef typename Kernel::Offset                           Offset;

  typedef typename Kernel::Construct_point_3                Construct_point_3;
  typedef typename Kernel::Construct_weighted_point_3       Construct_weighted_point_3;

public:
  typedef typename Functor::result_type                     result_type;

  Functor_with_offset_weighted_points_adaptor_3(const Functor_& functor,
                                                const Construct_point_3& cp,
                                                const Construct_weighted_point_3& wp)
    : Base(functor, cp), wp(wp)
  { }

  // gives access to calls with Point_3 arguments and without offset
  using Base::operator();

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Offset& o0, const Offset& o1) const {
    return operator()(wp(p0, o0), wp(p1, o1));
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2) const {
    return operator()(wp(p0, o0), wp(p1, o1), wp(p2, o2));
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2, const Offset& o3) const {
    return operator()(wp(p0, o0), wp(p1, o1), wp(p2, o2), wp(p3, o3));
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const Weighted_point_3& p4,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2, const Offset& o3,
                          const Offset& o4) const {
    return operator()(wp(p0, o0), wp(p1, o1), wp(p2, o2), wp(p3, o3), wp(p4, o4));
  }

  // for `Compare_power_distance_3`
  result_type operator() (const Point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2) const {
    return operator()(this->cp(p0, o0), wp(p1, o1), wp(p2, o2));
  }

  // for `Compare_weighted_squared_radius_3`
  result_type operator() (const Weighted_point_3& p0, const Offset& o0,
                          const FT w) const {
    return operator()(wp(p0, o0), w);
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Offset& o0, const Offset& o1,
                          const FT w) const {
    return operator()(wp(p0, o0), wp(p1, o1), w);
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2,
                          const FT w) const {
    return operator()(wp(p0, o0), wp(p1, o1), wp(p2, o2), w);
  }

  result_type operator() (const Weighted_point_3& p0, const Weighted_point_3& p1,
                          const Weighted_point_3& p2, const Weighted_point_3& p3,
                          const Offset& o0, const Offset& o1,
                          const Offset& o2, const Offset& o3,
                          const FT w) const {
    return operator()(wp(p0, o0), wp(p1, o1), wp(p2, o2), wp(p3, o3), w);
  }

  // for robust circumcenter_3
  result_type operator()(const Weighted_point_3& p0, const Weighted_point_3& p1,
                         const Weighted_point_3& p2, const Weighted_point_3& p3,
                         const Offset& o0, const Offset& o1,
                         const Offset& o2, const Offset& o3,
                         bool b) const {
    return operator()(wp(p0, o0), wp(p1, o1), wp(p2, o2), wp(p3, o3), b);
  }

  const Construct_weighted_point_3 wp;
};

}  // namespace CGAL

#endif /* CGAL_FUNCTOR_WITH_OFFSET_WEIGHTED_POINTS_ADAPTOR_3_H */
