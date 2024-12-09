// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Philipp Moeller
//

#ifndef CGAL_IS_RAY_INTERSECTION_GEOMTRAITS_H
#define CGAL_IS_RAY_INTERSECTION_GEOMTRAITS_H

#include <CGAL/license/AABB_tree.h>


#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/and.hpp>

namespace CGAL {
namespace internal {
namespace AABB_tree {

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_ray_3, Ray_3, false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_construct_source_3, Construct_source_3, false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_vector_3, Construct_vector_3, false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_construct_cartesian_const_iterator_3, Construct_cartesian_const_iterator_3, false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_cartesian_const_iterator_3, Cartesian_const_iterator_3, false)

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_ray_2, Ray_2, false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_construct_source_2, Construct_source_2, false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_vector_2, Construct_vector_2, false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_construct_cartesian_const_iterator_2, Construct_cartesian_const_iterator_2, false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_cartesian_const_iterator_2, Cartesian_const_iterator_2, false)

template<typename GeomTraits>
struct Is_ray_intersection_geomtraits_2
: std::bool_constant< Has_ray_2<GeomTraits>::value &&
                      Has_construct_source_2<GeomTraits>::value &&
                      Has_vector_2<GeomTraits>::value &&
                      Has_construct_cartesian_const_iterator_2<GeomTraits>::value &&
                      Has_cartesian_const_iterator_2<GeomTraits>::value >
{};

template<typename GeomTraits>
struct Is_ray_intersection_geomtraits
: std::bool_constant< Has_ray_3<GeomTraits>::value&&
                      Has_construct_source_3<GeomTraits>::value&&
                      Has_vector_3<GeomTraits>::value&&
                      Has_construct_cartesian_const_iterator_3<GeomTraits>::value&&
                      Has_cartesian_const_iterator_3<GeomTraits>::value >
{};


} // AABB_tree
} // internal
} // CGAL

#endif /* CGAL_IS_RAY_INTERSECTION_GEOMTRAITS_H */
