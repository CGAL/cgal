// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_MCFSKEL_UTILITY_H
#define CGAL_MCFSKEL_UTILITY_H

/// @cond CGAL_DOCUMENT_INTERNAL

/**
 * @file Utility.h
 * @brief This file contains some helper functions like splitting an edge at a
 * given point.
 */

#include <CGAL/Polygon_mesh_processing/measure.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <cmath>

namespace CGAL {
namespace internal {

template<class TriangleMesh, class TriangleMeshPointPMap, class Traits>
double get_surface_area(TriangleMesh& hg,
                        TriangleMeshPointPMap& hg_point_pmap,
                        const Traits& traits)
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  return PMP::area(hg
    , PMP::parameters::vertex_point_map(hg_point_pmap)
    .geom_traits(traits));
}

} //namespace internal
} //namespace CGAL

/// @endcond

#endif //CGAL_MCFSKEL_UTILITY_H
