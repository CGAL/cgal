// Copyright (c) 2018  Carnegie Mellon University (USA), GeometryFactory (France)
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s) : Christina Vaz, Keenan Crane, Andreas Fabri


#ifndef CGAL_INTRINSIC_DELAUNAY_TRIANGULATION_3_H
#define CGAL_INTRINSIC_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/license/Heat_method_3.h>

#include <CGAL/disable_warnings.h>
#include <set>

#include <CGAL/property_map.h>
#include <CGAL/double.h>
#include <Eigen/Cholesky>
#include <Eigen/Sparse>

#include <boost/foreach.hpp>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Dynamic_property_map.h>
#include <vector>
#include <CGAL/squared_distance_3.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intrinsic_Delaunay_Triangulation_3 {


    /**
     * Class `Intrinsic_Delaunay_Triangulation_3` is a ...
     * \tparam TriangleMesh a triangulated surface mesh, model of `FaceGraph` and `HalfedgeListGraph`
     * \tparam Traits a model of IntrinsicDelaunayTriangulation_3
     * \tparam VertexPointMap a model of `ReadablePropertyMap` with
     *        `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and
     *        `Traits::Point_3` as value type.
     *        The default is `typename boost::property_map< TriangleMesh, vertex_point_t>::%type`.
     *
     */
     template <typename TriangleMesh,
               typename Traits,
               typename EdgeLengthMap,
               typename FaceAreaMap,
               typename VertexPointMap = typename boost::property_map< TriangleMesh, vertex_point_t>::const_type,
               typename FaceIndexMap = typename boost::property_map< TriangleMesh, face_index_t>::const_type,
     class Intrinsic_Delaunay_Triangulation_3
     {

     };


} // namespace Intrinsic_Delaunay_Triangulation_3
} // namespace CGAL
#include <CGAL/enable_warnings.h>
#endif CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H
