// Copyright (c) 2024 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Stéphane Tayeb, Pierre Alliez, Camille Wormser
//

#ifndef CGAL_AABB_FILTERED_TRAITS_3_H
#define CGAL_AABB_FILTERED_TRAITS_3_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/AABB_traits_3.h>
#include <CGAL/Filtered_predicate.h>

namespace CGAL {

template<typename GeomTraits, typename AABBPrimitive, typename BboxMap = Default>
class AABB_filtered_traits_base_3 : public AABB_traits_base_3<GeomTraits, AABBPrimitive, BboxMap> {
  using Base = AABB_traits_base_3<GeomTraits, AABBPrimitive, BboxMap>;
public:
  typedef GeomTraits                                              Kernel;

  typedef typename Kernel::Exact_kernel                           EKernel;
  typedef typename Kernel::Approximate_kernel                     AKernel;
  typedef typename Kernel::C2E                                    C2E;
  typedef typename Kernel::C2F                                    C2F;

  // Exact traits is based on the exact kernel.
  typedef AABB_traits_base_3<EKernel, AABBPrimitive, BboxMap> Exact_traits;
  // Filtering traits is based on the filtering kernel.
  typedef AABB_traits_base_3<AKernel, AABBPrimitive, BboxMap> Filtering_traits;

  typedef Filtered_predicate<
    typename Exact_traits::Compare_distance,
    typename Filtering_traits::Compare_distance,
    C2E, C2F> Compare_distance;

  AABB_filtered_traits_base_3() : Base() {}
  AABB_filtered_traits_base_3(BboxMap bbm) : Base(bbm) {}

  Compare_distance compare_distance_object() const
  {
    typename Exact_traits::Compare_distance pe = Exact_traits().compare_distance_object();
    typename Filtering_traits::Compare_distance pf = Filtering_traits().compare_distance_object();

    return Compare_distance(pe, pf);
  }
};
}

#include <CGAL/AABB_tree/internal/AABB_statically_filtered_traits_3.h>

namespace CGAL {

template<typename GeomTraits, typename AABBPrimitive, typename BboxMap = Default, bool Has_static_filters = internal::Has_static_filters<GeomTraits>::value>
class AABB_filtered_traits_3;

template<typename GeomTraits, typename AABBPrimitive, typename BboxMap>
class AABB_filtered_traits_3<GeomTraits, AABBPrimitive, BboxMap, false> : public AABB_filtered_traits_base_3<GeomTraits, AABBPrimitive, BboxMap> {
  using Base = AABB_filtered_traits_base_3<GeomTraits, AABBPrimitive, BboxMap>;

public:
  AABB_filtered_traits_3() : Base() {}
  AABB_filtered_traits_3(BboxMap bbm) : Base(bbm) {}
};

template<typename GeomTraits, typename AABBPrimitive, typename BboxMap>
class AABB_filtered_traits_3<GeomTraits, AABBPrimitive, BboxMap, true> : public AABB_statically_filtered_traits_3<GeomTraits, AABBPrimitive, BboxMap> {
  using Base = AABB_statically_filtered_traits_3<GeomTraits, AABBPrimitive, BboxMap>;

public:
  AABB_filtered_traits_3() : Base() {}
  AABB_filtered_traits_3(BboxMap bbm) : Base(bbm) {}
};

} //namespace CGAL

#endif
