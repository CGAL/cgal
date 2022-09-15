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

#ifndef CGAL_MARCHING_CUBES_3_INTERNAL_MARCHING_CUBES_3_H
#define CGAL_MARCHING_CUBES_3_INTERNAL_MARCHING_CUBES_3_H

#include <CGAL/Isosurfacing_3/internal/Tables.h>
#include <tbb/concurrent_vector.h>

#include <array>
#include <atomic>
#include <bitset>
#include <cassert>
#include <mutex>

namespace CGAL {
namespace Isosurfacing {
namespace internal {

template <class Point_3, typename FT>
Point_3 vertex_interpolation(const Point_3& p0, const Point_3& p1, const FT d0, const FT d1, const FT iso_value) {

    FT mu;

    // don't divide by 0
    if (abs(d1 - d0) < 0.000001) {
        mu = 0.5;  // if both points have the same value, assume isolevel is in the middle
    } else {
        mu = (iso_value - d0) / (d1 - d0);
    }

    assert(mu >= 0.0 || mu <= 1.0);

    // linear interpolation
    return Point_3(p1.x() * mu + p0.x() * (1 - mu), p1.y() * mu + p0.y() * (1 - mu), p1.z() * mu + p0.z() * (1 - mu));
}

template <class Domain_, typename Corners_, typename Values_>
std::size_t get_cell_corners(const Domain_& domain, const typename Domain_::Cell_handle& cell,
                             const typename Domain_::FT iso_value, Corners_& corners, Values_& values) {

    typedef typename Domain_::Vertex_handle Vertex_handle;

    // collect function values and build index
    std::size_t v_id = 0;
    std::bitset<Domain_::VERTICES_PER_CELL> index = 0;
    for (const Vertex_handle& v : domain.cell_vertices(cell)) {
        // collect scalar values and computex index
        corners[v_id] = domain.position(v);
        values[v_id] = domain.value(v);

        if (values[v_id] >= iso_value) {
            index.set(v_id);
        }
        // next cell vertex
        v_id++;
    }

    return static_cast<std::size_t>(index.to_ullong());
}

template <class CellEdges, typename FT, typename Corners_, typename Values_, typename Vertices_>
void mc_construct_vertices(const CellEdges& cell_edges, const FT iso_value, const std::size_t i_case,
                           const Corners_& corners, const Values_& values, Vertices_& vertices) {

    // compute for this case the vertices
    std::size_t flag = 1;
    std::size_t e_id = 0;

    for (const auto& edge : cell_edges) {
        (void)edge;  // TODO

        if (flag & Cube_table::intersected_edges[i_case]) {

            // generate vertex here, do not care at this point if vertex already exist

            const int v0 = Cube_table::edge_to_vertex[e_id][0];
            const int v1 = Cube_table::edge_to_vertex[e_id][1];

            vertices[e_id] = vertex_interpolation(corners[v0], corners[v1], values[v0], values[v1], iso_value);
        }
        flag <<= 1;
        e_id++;
    }
}

template <typename Vertices_, class TriangleList>
void mc_construct_triangles(const int i_case, const Vertices_& vertices, TriangleList& triangles) {
    // construct triangles
    for (int t = 0; t < 16; t += 3) {

        const int t_index = i_case * 16 + t;
        // if (e_tris_list[t_index] == 0x7f)
        if (Cube_table::triangle_cases[t_index] == -1) break;

        const int eg0 = Cube_table::triangle_cases[t_index + 0];  // TODO: move more of this stuff into the table
        const int eg1 = Cube_table::triangle_cases[t_index + 1];
        const int eg2 = Cube_table::triangle_cases[t_index + 2];

        // insert new triangle in list
        std::array<Point, 3> points;
        points[0] = vertices[eg0];
        points[1] = vertices[eg1];
        points[2] = vertices[eg2];

        triangles.push_back(points);
    }
}

template <class TriangleList, class PointRange, class PolygonRange>
void to_indexed_face_set(const TriangleList& triangle_list, PointRange& points, PolygonRange& polygons) {
    for (auto& triangle : triangle_list) {
        const std::size_t id = points.size();

        points.push_back(triangle[0]);
        points.push_back(triangle[1]);
        points.push_back(triangle[2]);

        polygons.push_back({});
        auto& triangle = polygons.back();
        triangle.push_back(id + 2);
        triangle.push_back(id + 1);
        triangle.push_back(id + 0);
    }
}

template <class Domain_>
class Marching_cubes_functor {
private:
    typedef Domain_ Domain;
    typedef typename Domain::FT FT;
    typedef typename Domain::Point Point;
    typedef typename Domain::Cell_handle Cell_handle;

    typedef tbb::concurrent_vector<std::array<Point, 3>> Triangle_list;

public:
    Marching_cubes_functor(const Domain& domain, const FT iso_value) : domain(domain), iso_value(iso_value) {}


    void operator()(const Cell_handle& cell) {

        assert(domain.cell_vertices(cell).size() == 8);
        assert(domain.cell_edges(cell).size() == 12);

        FT values[8];
        Point corners[8];
        const int i_case = get_cell_corners(domain, cell, iso_value, corners, values);

        const int all_bits_set = (1 << (8 + 1)) - 1;  // last 8 bits are 1
        if (i_case == 0 || i_case == all_bits_set) {
            return;
        }

        std::array<Point, 12> vertices;
        mc_construct_vertices(domain.cell_edges(cell), iso_value, i_case, corners, values, vertices);

        mc_construct_triangles(i_case, vertices, triangle_list);
    }

    const Triangle_list& get_triangles() const {
        return triangle_list;
    }

private:
    const Domain& domain;
    FT iso_value;

    Triangle_list triangle_list;
};

}  // namespace internal
}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_MARCHING_CUBES_3_INTERNAL_MARCHING_CUBES_3_H
