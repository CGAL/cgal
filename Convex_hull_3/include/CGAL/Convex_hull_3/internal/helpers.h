// Copyright (c) 2025 Geometry Factory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Léo Valque
//

#ifndef CGAL_CONVEX_HULL_INTERNAL_HELPERS_H
#define CGAL_CONVEX_HULL_INTERNAL_HELPERS_H

#include <CGAL/IO/helpers.h>
#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/Container_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Convex_hull_hierarchy.h>
#include <CGAL/extreme_point_3.h>

#ifdef CGAL_PROFILE_CONVEX_HULL_DO_INTERSECT
      size_t nb_visited=0;
#endif

namespace CGAL::Convex_hull_3::internal{

template <class T>
inline constexpr bool is_instance_of_CHH = false;

template <class Arg>
inline constexpr bool is_instance_of_CHH< Convex_hull_hierarchy<Arg> > = true;

// template class to deduce the GT from Convex and NamedParameters
template<class Convex,
         class NamedParameters,
         bool Is_range=CGAL::IO::internal::is_Range_v<Convex> >
struct GetGeomTraitsFromConvex{
  typedef typename GetGeomTraits<Convex, NamedParameters>::type type;
};

template<class Convex, class NamedParameters>
struct GetGeomTraitsFromConvex<Convex, NamedParameters, true>{
  typedef typename Point_set_processing_3_np_helper<Convex, NamedParameters>::Geom_traits type;
};

template<class Mesh, class NamedParameters>
struct GetGeomTraitsFromConvex<Convex_hull_hierarchy<Mesh>, NamedParameters, false>{
  typedef typename GetGeomTraits<Mesh, NamedParameters>::type type;
};

template <class Convex, class Direction_3, class NamedParameters>
typename Kernel_traits<Direction_3>::Kernel::Point_3 extreme_point_3_wrapper(const Convex &c, const Direction_3 &dir, const NamedParameters &np){
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  if constexpr(Convex_hull_3::internal::is_instance_of_CHH<Convex>){
    return c.extreme_point_3(dir, np);
  } else if constexpr(CGAL::IO::internal::is_Range_v<Convex>){
    using NP_helper= Point_set_processing_3_np_helper<Convex, NamedParameters>;

    using Point_GT = typename NP_helper::Geom_traits;
    using Default_GT = typename Kernel_traits<Direction_3>::Kernel;
    using GT=typename internal_np::Lookup_named_param_def <
        internal_np::geom_traits_t,
        NamedParameters,
        Default_GT
      > ::type;

    using Default_geom_traits_converter = Cartesian_converter<Point_GT, GT>;
    using GTC=typename internal_np::Lookup_named_param_def <
        internal_np::geom_traits_converter_t,
        NamedParameters,
        Default_geom_traits_converter
      > ::type;
    GTC converter = choose_parameter<GTC>(get_parameter(np, internal_np::geom_traits_converter));

    auto point_map = NP_helper::get_const_point_map(c, np);
    return converter(get(point_map, extreme_point_3(c, dir, np)));
  } else { // Convex is a Graph
    using GetGeomTraits = GetGeomTraits<Convex, NamedParameters>;
    using GraphGT= typename GetGeomTraits::type;

    using Default_GT = typename Kernel_traits<Direction_3>::Kernel;
    using GT=typename internal_np::Lookup_named_param_def <
      internal_np::geom_traits_t,
      NamedParameters,
      Default_GT
    > ::type;

    using Default_geom_traits_converter = Cartesian_converter<GraphGT, GT>;
    using GTC=typename internal_np::Lookup_named_param_def <
        internal_np::geom_traits_converter_t,
        NamedParameters,
        Default_geom_traits_converter
      > ::type;
    GTC converter = choose_parameter<GTC>(get_parameter(np, internal_np::geom_traits_converter));
    auto point_map = GetVertexPointMap<Convex, NamedParameters>::get_const_map(np, c);
    return converter(get(point_map, extreme_vertex_3(c, dir, np)));
  }
}

}

#endif // CGAL_CONVEX_HULL_INTERNAL_HELPERS_H