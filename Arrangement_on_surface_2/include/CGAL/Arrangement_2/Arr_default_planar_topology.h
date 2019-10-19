// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Eric Berberich    <ericb@post.tau.ac.il>

#ifndef CGAL_ARR_DEFAULT_PLANAR_TOPOLOGY_H
#define CGAL_ARR_DEFAULT_PLANAR_TOPOLOGY_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of a helper template that uses partial specialization to
 * select the default planar topology-traits class best fit for a
 * given geometry-traits class.
 */

#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_bounded_planar_topology_traits_2.h>
#include <CGAL/Arr_unb_planar_topology_traits_2.h>

namespace CGAL {

namespace internal {

template <class GeomTraits, class Dcel, class Tag>
struct Default_planar_topology_impl {};

template <class GeomTraits, class Dcel>
struct Default_planar_topology_impl<GeomTraits, Dcel,
                                    Arr_all_sides_oblivious_tag>
{
  // A topology-traits class that supports only bounded curves:
  typedef Arr_bounded_planar_topology_traits_2<GeomTraits, Dcel>    Traits;

};

template <class GeomTraits, class Dcel>
struct Default_planar_topology_impl<GeomTraits, Dcel,
                                    Arr_not_all_sides_oblivious_tag>
{
  // A topology-traits class that supports unbounded curves:
  typedef Arr_unb_planar_topology_traits_2<GeomTraits, Dcel>        Traits;
};

} // namespace internal

template <typename GeomTraits, typename Dcel >
struct Default_planar_topology :
  public internal::Default_planar_topology_impl<
    GeomTraits, Dcel,
    typename Arr_are_all_sides_oblivious_tag<
      typename internal::Arr_complete_left_side_category<GeomTraits>::Category,
      typename internal::Arr_complete_bottom_side_category<GeomTraits>::Category,
      typename internal::Arr_complete_top_side_category<GeomTraits>::Category,
      typename internal::Arr_complete_right_side_category<GeomTraits>::Category>::result
>
{};

} //namespace CGAL

#endif
