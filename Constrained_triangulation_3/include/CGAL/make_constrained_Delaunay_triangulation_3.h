// Copyright (c) 2024  GeometryFactory Sarl (France).
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
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MAKE_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
#define CGAL_MAKE_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/license/Constrained_triangulation_3.h>

#include <CGAL/Constrained_Delaunay_triangulation_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_cell_data_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_vertex_base_3.h>

namespace CGAL {

/*!
 * \ingroup PkgCT_3Functions
 * \brief Create a 3D constrained Delaunay triangulation conforming to the faces of a polygon mesh.
 *
 * The polygon mesh represents the polygonal constraints that will be enforced during the triangulation process.
 *
 * By default, each face of the polygon mesh is considered as a polygonal constraint for the triangulation. The
 * named parameter `face_patch_map` can be used to describe bigger polygonal constraints, possibly with holes. If
 * used, the argument of that parameter must be a property map that maps each face of the polygon mesh to a patch
 * identifier. Faces with the same patch identifier are considered as part of the same surface patch. Each of those
 * surface patches (defined as the union of the mesh faces with a given patch id) is supposed to be a polygon or a
 * polygon with holes, with coplanar vertices (or almost coplanar up to the precision of the number type used).
 *
 * The generated triangulation will be constrained to conform to the faces of the polygon mesh, or to the surface patches
 * described by the `face_patch_map` property map if provided.
 *
 * \tparam Triangulation_3 An instance of the `Triangulation_3` class template.
 *         - Its `Geom_traits` type must be a model of `ConstrainedDelaunayTriangulationTraits_3`,
 *         - Its point type must be constructible from the point type of the polygon mesh,
 *         - its `Vertex` type must be a model of `ConstrainedDelaunayTriangulationVertexBase_3`, and
 *         - its `Cell` type must be a model of `ConstrainedDelaunayTriangulationCellBase_3`.
 * \tparam PolygonMesh a model of `FaceListGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param mesh The polygon mesh representing the constraints.
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `mesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Traits::Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `PolygonMesh`.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{`Traits`}
 *     \cgalParamDefault{the default constructed traits object `Traits{}`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{face_patch_map}
 *    \cgalParamDescription{a property map associating a patch identifier to each face of `mesh`}
 *    \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
 *                   as key type and with any value type that is a *regular* type}
 *   \cgalParamExtra{If this parameter is omitted, each face of the mesh is considered as a separate patch.}
 *   \cgalParamExtra{Faces with the same patch identifier are considered as part of the same surface patch.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre `mesh` must not have self-intersections:
 *      \link CGAL::Polygon_mesh_processing::does_self_intersect
 *      `CGAL::Polygon_mesh_processing::does_self_intersect(mesh, np) == false`
 *      \endlink
 *
 * \todo Create a documentation page to describe the concept *regular*, and link it to https://en.cppreference.com/w/cpp/concepts/regular
 */
template <typename Triangulation_3, typename PolygonMesh, typename NamedParams = parameters::Default_named_parameters>
Triangulation_3
make_constrained_Delaunay_triangulation_3(const PolygonMesh& mesh, const NamedParams& np = parameters::default_values())
{
  Constrained_Delaunay_triangulation_3<typename Triangulation_3::Geom_traits, Triangulation_3> cdt(mesh, np);
  Triangulation_3 tr = std::move(cdt).triangulation();

  for(auto vh : tr.all_vertex_handles()) {
    vh->sync();
  }

  for(auto ch : tr.all_cell_handles()) {
    ch->set_subdomain_index(1);
  }

  std::stack<typename Triangulation_3::Cell_handle> stack;
  stack.push(tr.infinite_cell());
  while(!stack.empty()) {
    auto ch = stack.top();
    stack.pop();
    ch->set_subdomain_index(0);
    for(int i = 0; i < 4; ++i) {
      if(ch->is_facet_on_surface(i)) continue;
      auto n = ch->neighbor(i);
      if(n->subdomain_index() == 1) {
        stack.push(n);
      }
    }
  }

  return tr;
}

} // end namespace CGAL

#endif // CGAL_MAKE_CONSTRAINED_DELAUNAY_TRIANGULATION_3_H
