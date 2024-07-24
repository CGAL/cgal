// Copyright (c) 2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: ( GPL-3.0-or-later OR LicenseRef-Commercial ) AND MIT
//
// Author(s)     : Julian Stahl
//
// This file incorporates work covered by the following copyright and permission notice:
//
//     MIT License
//
//     Copyright (c) 2020 Roberto Grosso
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy
//     of this software and associated documentation files (the "Software"), to deal
//     in the Software without restriction, including without limitation the rights
//     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//     copies of the Software, and to permit persons to whom the Software is
//     furnished to do so, subject to the following conditions:
//
//     The above copyright notice and this permission notice shall be included in all
//     copies or substantial portions of the Software.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//     SOFTWARE.
//
//
// The code below uses the version of
// https://github.com/rogrosso/tmc available on 15th of September 2022.
//

#ifndef CGAL_ISOSURFACING_3_INTERNAL_MARCHING_CUBES_FUNCTORS_H
#define CGAL_ISOSURFACING_3_INTERNAL_MARCHING_CUBES_FUNCTORS_H

#include <CGAL/license/Isosurfacing_3.h>

#include <CGAL/Isosurfacing_3/internal/tables.h>

#include <CGAL/assertions.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/enumerable_thread_specific.h>
#else
#include <vector>
#endif

#include <array>
#include <bitset>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

// Interpolate linearly between two vertex locations v0, v1 with values d0 and d1 according to the isovalue
template <typename GeomTraits>
typename GeomTraits::Point_3 vertex_interpolation(const typename GeomTraits::Point_3& p0,
                                                  const typename GeomTraits::Point_3& p1,
                                                  const typename GeomTraits::FT d0,
                                                  const typename GeomTraits::FT d1,
                                                  const typename GeomTraits::FT isovalue,
                                                  const GeomTraits& gt)
{
  using FT = typename GeomTraits::FT;

  typename GeomTraits::Compute_x_3 x_coord = gt.compute_x_3_object();
  typename GeomTraits::Compute_y_3 y_coord = gt.compute_y_3_object();
  typename GeomTraits::Compute_z_3 z_coord = gt.compute_z_3_object();
  typename GeomTraits::Construct_point_3 point = gt.construct_point_3_object();

  FT mu = FT(0);

  // @todo, technically we should be using the edge intersection oracle here, but there is a nuance
  // between MC and DC on the handling of edges that have val0 = val1 = isovalue: in MC we assume
  // the isosurface is in the middle, in DC we assume the isosurface is not intersecting the edge.
  // In the oracle, we follow DC right now. Could put a Boolean parameter, but it's ugly.

  // don't divide by 0
  if(abs(d1 - d0) < 0.000001) // @fixme hardcoded bound
    mu = FT(0.5);  // if both points have the same value, assume isolevel is in the middle
  else
    mu = (isovalue - d0) / (d1 - d0);

  CGAL_assertion(mu >= FT(0.0) || mu <= FT(1.0));

  // linear interpolation
  return point(x_coord(p1) * mu + x_coord(p0) * (FT(1) - mu),
               y_coord(p1) * mu + y_coord(p0) * (FT(1) - mu),
               z_coord(p1) * mu + z_coord(p0) * (FT(1) - mu));
}

// retrieves the values of a cell and return the lookup index
// if the cell is completely above or below the isovalue, corner points are not computed
template <typename Domain,
          typename Corners,
          typename Values>
std::size_t get_cell_corners(const Domain& domain,
                             const typename Domain::cell_descriptor& cell,
                             const typename Domain::Geom_traits::FT isovalue,
                             Corners& corners,
                             Values& values)
{
  using vertex_descriptor = typename Domain::vertex_descriptor;

  const auto& vertices = domain.cell_vertices(cell);

  // collect function values and build index
  std::size_t v_id = 0;
  std::bitset<Domain::VERTICES_PER_CELL> index = 0;
  for(const vertex_descriptor& v : vertices)
  {
    values[v_id] = domain.value(v);
    if(values[v_id] >= isovalue)
      index.set(v_id);

    ++v_id;
  }

  if(index.all() || index.none()) // nothing's happening in this cell
    return static_cast<std::size_t>(index.to_ullong());

  v_id = 0;
  for(const vertex_descriptor& v : vertices)
    corners[v_id++] = domain.point(v);

  return static_cast<std::size_t>(index.to_ullong());
}

// creates the vertices on the edges of one cell
template <typename Corners,
          typename Values,
          typename Domain,
          typename Vertices>
void MC_construct_vertices(const typename Domain::cell_descriptor& cell,
                           const std::size_t i_case,
                           const Corners& corners,
                           const Values& values,
                           const typename Domain::Geom_traits::FT isovalue,
                           const Domain& domain,
                           Vertices& vertices)
{
  using Cell_edges = typename Domain::Cell_edges;
  using edge_descriptor = typename Domain::edge_descriptor;

  const Cell_edges& cell_edges = domain.cell_edges(cell);

  // compute for this case the vertices
  std::size_t flag = 1;
  std::size_t e_id = 0;

  for(const edge_descriptor& e : cell_edges)
  {
    CGAL_USE(e);

    if(flag & Cube_table::intersected_edges[i_case])
    {
      // generate vertex here, do not care at this point if vertex already exists
      const int v0 = Cube_table::edge_to_vertex[e_id][0];
      const int v1 = Cube_table::edge_to_vertex[e_id][1];

      vertices[e_id] = vertex_interpolation(corners[v0], corners[v1],
                                            values[v0], values[v1],
                                            isovalue, domain.geom_traits());
    }

    flag <<= 1;
    ++e_id;
  }
}

// connects the vertices of one cell to form triangles
template <typename Vertices,
          typename TriangleList>
void MC_construct_triangles(const std::size_t i_case,
                            const Vertices& vertices,
                            TriangleList& triangles)
{
  // construct triangles
  for(int t=0; t<16; t+=3)
  {
    const std::size_t t_index = i_case * 16 + t;

    // if(e_tris_list[t_index] == 0x7f)
    if(Cube_table::triangle_cases[t_index] == -1)
      break;

    // @todo move more of this stuff into the table
    const int eg0 = Cube_table::triangle_cases[t_index + 0];
    const int eg1 = Cube_table::triangle_cases[t_index + 1];
    const int eg2 = Cube_table::triangle_cases[t_index + 2];

    // insert new triangle in list
#ifdef CGAL_LINKED_WITH_TBB
    auto& tris = triangles.local();
#else
    auto& tris = triangles;
#endif

    tris.push_back({vertices[eg0], vertices[eg1], vertices[eg2]});
  }
}

template <typename TriangleRange,
          typename PointRange,
          typename PolygonRange>
void triangles_to_polygon_soup(const TriangleRange& triangles,
                               PointRange& points,
                               PolygonRange& polygons)
{
#ifdef CGAL_LINKED_WITH_TBB
    for(const auto& triangle_list : triangles)
    {
#else
      const auto& triangle_list = triangles;
#endif

      for(const auto& triangle : triangle_list)
      {
        const std::size_t id = points.size();

        points.push_back(triangle[0]);
        points.push_back(triangle[1]);
        points.push_back(triangle[2]);

        // simply use increasing indices
        polygons.push_back({id + 2, id + 1, id + 0});

        // just a safeguard against arrays of the wrong size
        CGAL_assertion(polygons.back().size() == 3);
      }

#ifdef CGAL_LINKED_WITH_TBB
    }
#endif
}

// Marching Cubes implemented as a functor that runs on every cell of the grid
template <typename Domain_>
class Marching_cubes_3
{
public:
  using Domain = Domain_;

  using Geom_traits = typename Domain::Geom_traits;
  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;

  using cell_descriptor = typename Domain::cell_descriptor;

#ifdef CGAL_LINKED_WITH_TBB
  using Triangles = tbb::enumerable_thread_specific<std::vector<std::array<Point_3, 3>>>;
#else
  using Triangles = std::vector<std::array<Point_3, 3> >;
#endif

private:
  const Domain& m_domain;
  const FT m_isovalue;

  Triangles m_triangles;

public:
  // creates a Marching Cubes functor for a domain and isovalue
  Marching_cubes_3(const Domain& domain,
                   const FT isovalue)
    : m_domain(domain),
      m_isovalue(isovalue)
  { }

  // returns the created triangle list
  Triangles& triangles()
  {
    return m_triangles;
  }

public:
  // computes one cell
  void operator()(const cell_descriptor& cell)
  {
    CGAL_precondition(m_domain.cell_vertices(cell).size() == 8);
    CGAL_precondition(m_domain.cell_edges(cell).size() == 12);

    // @todo for SDFs, we could query at the center of the voxel an early exit

    constexpr std::size_t vpc = Domain::VERTICES_PER_CELL;

    std::array<FT, vpc> values;
    std::array<Point_3, vpc> corners;
    const std::size_t i_case = get_cell_corners(m_domain, cell, m_isovalue, corners, values);

    // skip empty / full cells
    constexpr std::size_t ones = (1 << vpc) - 1;
    if((i_case & ones) == ones || // all bits set
       (i_case & ones) == 0) // no bits set
      return;

    std::array<Point_3, 12> vertices;
    MC_construct_vertices(cell, i_case, corners, values, m_isovalue, m_domain, vertices);

    MC_construct_triangles(i_case, vertices, m_triangles);
  }
};

} // namespace internal
} // namespace Isosurfacing
} // namespace CGAL

#endif // CGAL_ISOSURFACING_3_INTERNAL_MARCHING_CUBES_FUNCTORS_H
