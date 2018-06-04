// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_HAS_MEMBER_VISITED_H
#define CGAL_HAS_MEMBER_VISITED_H

#include <boost/mpl/has_xxx.hpp>

namespace CGAL {
namespace internal {

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_member_visited,
                                    Has_visited_for_vertex_extractor,
                                    false)

} // end namespace internal
} // end namespace CGAL

#endif // CGAL_HAS_NESTED_TYPE_BARE_POINT_H
