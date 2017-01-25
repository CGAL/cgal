// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_ray_3,Ray_3,false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_construct_source_3,Construct_source_3,false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_vector_3,Construct_vector_3,false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_construct_cartesian_const_iterator_3,Construct_cartesian_const_iterator_3,false)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_cartesian_const_iterator_3,Cartesian_const_iterator_3,false)

template<typename GeomTraits>
struct Is_ray_intersection_geomtraits
: boost::mpl::and_< Has_ray_3<GeomTraits>,
                    Has_construct_source_3<GeomTraits>,
                    Has_vector_3<GeomTraits>,
                    Has_construct_cartesian_const_iterator_3<GeomTraits>,
                    Has_cartesian_const_iterator_3<GeomTraits> >::type
{};


} // AABB_tree
} // internal
} // CGAL

#endif /* CGAL_IS_RAY_INTERSECTION_GEOMTRAITS_H */
