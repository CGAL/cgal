// Copyright (c) 2021 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Antonio Gomes, Dmitry Anisimov
//

#ifndef CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_3_H
#define CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_3_H

// #include <CGAL/license/Barycentric_coordinates_3.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>

namespace CGAL {
namespace Barycentric_coordinates {

  template<
  typename PolygonMesh,
  typename GeomTraits,
  typename VertexToPointMap = typename property_map_selector<PolygonMesh, CGAL::vertex_point_t>::const_type>
  class Wachspress_coordinates_3 {

  public:
    using Polygon_mesh = PolygonMesh;
    using Geom_traits = GeomTraits;
    using Vertex_to_point_map = VertexToPointMap;

		typedef typename GeomTraits::FT FT;
    typedef typename GeomTraits::Point_3 Point_3;

  public:
    Wachspress_coordinates_3(
    	const PolygonMesh& polygon_mesh,
    	const Computation_policy_3 policy,
    	const VertexToPointMap vertex_to_point_map,
    	const GeomTraits traits = GeomTraits()) :
    m_polygon_mesh(polygon_mesh),
    m_computation_policy(policy),
    m_vertex_to_point_map(vertex_to_point_map),
    m_traits(traits) {

    	// preconditions, resize containers, etc.

    	m_weights.resize(vertices(m_polygon_mesh).size());
    }

    Wachspress_coordinates_3(
    	const PolygonMesh& polygon_mesh,
    	const Computation_policy_3 policy =
    	Computation_policy_3::DEFAULT,
    	const GeomTraits traits = GeomTraits()) :
    Wachspress_coordinates_3(
      polygon_mesh,
      policy,
      get_const_property_map(CGAL::vertex_point, polygon_mesh),
      traits) { }

    template<typename OutIterator>
    OutIterator operator()(const Point_3& query, OutIterator c_begin) {
      return compute(query, c_begin);
    }

  private:
  	const PolygonMesh& m_polygon_mesh;
  	const Computation_policy_3 m_computation_policy; // skip it for the moment
  	const VertexToPointMap m_vertex_to_point_map; // use it to map vertex to Point_3
  	const GeomTraits m_traits;

  	std::vector<FT> m_weights;

  	// put here any other necessary global containers

  	template<typename OutputIterator>
    OutputIterator compute(
      const Point_3& query, OutputIterator coordinates) {

    	// Compute weights.
    	const FT sum = compute_weights(query);
    	CGAL_assertion(sum > FT(0));

      // The coordinates must be saved in the same order as vertices in the vertex range.
      std::size_t vi = 0;
      const auto vd = vertices(m_polygon_mesh);
      for (const auto vertex : vd) {
      	CGAL_assertion(vi < m_weights.size());
      	const FT coordinate = m_weights[vi] / sum;
        *(coordinates++) = coordinate;
        ++vi;
      }
    }

    FT compute_weights(const Point_3& query) {

    	// Sum of weights to normalize them later.
    	FT sum = FT(0);
			// Vertex index.
    	std::size_t vi = 0;
    	// Vertex range, you can make it global.
    	const auto vd = vertices(m_polygon_mesh);
    	for (const auto vertex : vd) {

            std::cout << get(m_vertex_to_point_map, vertex) << "\n";
    		const FT weight = FT(vi); // compute it here for query
    		CGAL_assertion(vi < m_weights.size());
    		m_weights[vi] = weight;
    		sum += weight;
    		++vi; // update vi
      }
      CGAL_assertion(sum > FT(0));
      return sum;
    }
  };

  template<
  typename Point_3,
  typename OutIterator,
  typename GeomTraits>
  OutIterator wachspress_coordinates_3(
    const CGAL::Surface_mesh<Point_3>& surface_mesh,
    const Point_3& query,
    OutIterator c_begin,
    const Computation_policy_3 policy =
    Computation_policy_3::DEFAULT) {

    using Geom_Traits = typename Kernel_traits<Point_3>::Kernel;
    using SM = CGAL::Surface_mesh<Point_3>;

    Wachspress_coordinates_3<SM, Geom_Traits> wachspress(surface_mesh, policy);
    return wachspress(query, c_begin);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_3_H
