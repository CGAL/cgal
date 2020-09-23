// Copyright (c) 2005-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                  Laurent Rineau <Laurent.Rineau@geometryfactory.com>

// This traits override Construct_circumcenter_3 and
// Construct_weighted_circumcenter_3
// to get robust ones when called on slivers

#ifndef CGAL_ROBUST_CIRCUMCENTER_TRAITS_3_H
#define CGAL_ROBUST_CIRCUMCENTER_TRAITS_3_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/number_utils_classes.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Robust_construction.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace CGAL {

template < typename K >
class Robust_construct_weighted_circumcenter_3
{
public:
  typedef Exact_predicates_exact_constructions_kernel EK;

  typedef typename K::Weighted_point_3               Weighted_point_3;
  typedef typename K::Bare_point                     Bare_point;
  typedef typename K::FT                             FT;

  typedef Bare_point       result_type;

  typedef Cartesian_converter<typename K::Kernel, EK>  To_exact;
  typedef Cartesian_converter<EK, typename K::Kernel>  Back_from_exact;


  Bare_point operator() ( const Weighted_point_3 & p,
                          const Weighted_point_3 & q,
                          const Weighted_point_3 & r,
                          const Weighted_point_3 & s) const
  {
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Construct_weighted_circumcenter_3
      exact_weighted_circumcenter = EK().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q),
                                                       to_exact(r),
                                                       to_exact(s)));
  }

  Bare_point operator() ( const Weighted_point_3 & p,
                          const Weighted_point_3 & q,
                          const Weighted_point_3 & r) const
  {
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Construct_weighted_circumcenter_3
      exact_weighted_circumcenter = EK().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q),
                                                       to_exact(r)));
  }

  Bare_point operator() ( const Weighted_point_3 & p,
                          const Weighted_point_3 & q) const
  {
    To_exact to_exact;
    Back_from_exact Back_from_exact;
    EK::Construct_weighted_circumcenter_3
      exact_weighted_circumcenter = EK().construct_weighted_circumcenter_3_object();

    return back_from_exact(exact_weighted_circumcenter(to_exact(p),
                                                       to_exact(q)));
  }
};

template < class K>
class Robust_circumcenter_traits_3
  : public K
{
  typedef Exact_predicates_exact_constructions_kernel EK;
 public:
  typedef CGAL::Robust_construction<EK::Construct_circumcenter_3,
                                    Cartesian_converter<K, EK>,
                                    Cartesian_converter<EK, K>,
                                    typename K::Point_3 >   Construct_circumcenter_3;
  typedef CGAL::Robust_construction<EK::Compute_squared_radius_3,
                                    Cartesian_converter<K, EK>,
                                    Cartesian_converter<EK, K>,
                                    typename K::FT >        Compute_squared_radius_3;

  Construct_circumcenter_3
  construct_circumcenter_3_object() const
  { return Construct_circumcenter_3(); }

  Compute_squared_radius_3
  compute_squared_radius_3_object() const
  { return Compute_squared_radius_3(); }
};

template < class K>
class Robust_weighted_circumcenter_traits_3
  : public K
{
  typedef Exact_predicates_exact_constructions_kernel EK;
 public:
  typedef typename Robust_circumcenter_traits_3<typename K::Kernel>::Construct_circumcenter_3 Construct_circumcenter_3;
  typedef typename Robust_circumcenter_traits_3<K>::Compute_squared_radius_3                  Compute_squared_radius_3;

  Construct_circumcenter_3
  construct_circumcenter_3_object() const
  { return Construct_circumcenter_3(); }

  typedef CGAL::Robust_construct_weighted_circumcenter_3<K> Construct_weighted_circumcenter_3;

  Construct_weighted_circumcenter_3
  construct_weighted_circumcenter_3_object() const
  { return Construct_weighted_circumcenter_3(); }

  Compute_squared_radius_3
  compute_squared_radius_3_object() const
  { return Compute_squared_radius_3(); }
};

} //namespace CGAL

#endif //CGAL_ROBUST_CIRCUMCENTER_TRAITS_3_H
