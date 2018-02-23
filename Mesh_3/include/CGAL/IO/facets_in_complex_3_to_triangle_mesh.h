// Copyright (c) 2009-2017 GeometryFactory (France).
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
// Author(s)     : Maxime Gimeno,
//                 Mael Rouxel-Labbé

#ifndef CGAL_FACETS_IN_COMPLEX_3_TO_TRIANGLE_MESH_H
#define CGAL_FACETS_IN_COMPLEX_3_TO_TRIANGLE_MESH_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/array.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Hash_handles_with_or_without_timestamps.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>

#include <cstddef>
#include <iterator>
#include <vector>

namespace CGAL {

//! \ingroup PkgMesh_3Functions
//!
//! \brief builds a `TriangleMesh` from the surface facets, with a consistent orientation
//!        at the interface of two subdomains.
//!
//! This function exports the surface as a `TriangleMesh` and appends it to `graph`, using
//! `orient_polygon_soup()`.
//!
//! @tparam C3T3 a model of `MeshComplexWithFeatures_3InTriangulation_3`.
//! @tparam TriangleMesh a model of `MutableFaceGraph` with an internal point property map.
//!         The point type should be compatible with the one used in `C3T3`.
//!
//! @param c3t3 an instance of `C3T3`.
//! @param graph an instance of `TriangleMesh`.
template<class C3T3, class TriangleMesh>
void facets_in_complex_3_to_triangle_mesh(const C3T3& c3t3, TriangleMesh& graph)
{
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type  VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type              Point_3;

  typedef typename C3T3::Triangulation                                   Tr;

  typedef typename Tr::Vertex_handle                                     Vertex_handle;
  typedef typename Tr::Cell_handle                                       Cell_handle;
  typedef typename Tr::Weighted_point                                    Weighted_point;

  typedef typename C3T3::Facets_in_complex_iterator                      Ficit;

  typedef CGAL::Hash_handles_with_or_without_timestamps                  Hash_fct;
  typedef boost::unordered_map<Vertex_handle, std::size_t, Hash_fct>     VHmap;

  typedef CGAL::cpp11::array<std::size_t, 3>                             Face;

  typename std::iterator_traits<Ficit>::difference_type nf =
    std::distance(c3t3.facets_in_complex_begin(), c3t3.facets_in_complex_end());

  std::vector<Face> faces;
  faces.reserve(nf);
  std::vector<Point_3> points;
  points.reserve(nf / 2); // approximating Euler

  VHmap vh_to_ids;
  std::size_t inum = 0;

  for(Ficit fit = c3t3.facets_in_complex_begin(),
            end = c3t3.facets_in_complex_end(); fit != end; ++fit)
  {
    Cell_handle c = fit->first;
    int s = fit->second;
    Face f;

    for(std::size_t i=1; i<4; ++i)
    {
      typename VHmap::iterator map_entry;
      bool is_new;
      std::size_t id(-1);
      Vertex_handle v = c->vertex((s+i)&3);
      CGAL_assertion(v != Vertex_handle() && !c3t3.triangulation().is_infinite(v));

      boost::tie(map_entry, is_new) = vh_to_ids.insert(std::make_pair(v, inum));
      if(is_new)
      {
        const Weighted_point& p = c3t3.triangulation().point(c, (s+i)&3);
        points.push_back(Point_3(CGAL::to_double(p.x()),
                                 CGAL::to_double(p.y()),
                                 CGAL::to_double(p.z())));
        id = inum++;
      }
      else
      {
        id = map_entry->second;
      }

      f[i-1] = id;
    }

    faces.push_back(f);
  }

  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces);

  CGAL_assertion(CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(faces));
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, graph);
}

} // namespace CGAL
#endif // CGAL_FACETS_IN_COMPLEX_3_TO_TRIANGLE_MESH_H
