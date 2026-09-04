// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : François Protais

#ifndef CGAL_MESH_SMOOTHING_3_INTERNAL_HILBERT_SORT_H
#define CGAL_MESH_SMOOTHING_3_INTERNAL_HILBERT_SORT_H

#include <CGAL/license/Mesh_smoothing_3.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>

namespace CGAL {

namespace Mesh_smoothing_3_internal {

struct Mesh_permutation
{
    std::vector<std::size_t> vertex_new_to_old;
    std::vector<std::size_t> vertex_old_to_new;
    std::vector<std::size_t> tet_new_to_old;
    std::vector<std::size_t> tet_old_to_new;
};

// gpt/copilot produced
Mesh_permutation hilbert_sort_mesh(
    std::vector<double>& points,
    std::vector<std::array<unsigned, 4>>& tets)
{
    using Kernel  = CGAL::Simple_cartesian<double>;
    using Point_3 = Kernel::Point_3;

    Mesh_permutation permutation;

    // -------------------------------------------------------------------------
    // Vertices
    // -------------------------------------------------------------------------

    std::vector<Point_3> cgal_points;
    cgal_points.reserve(points.size() / 3);

    for (std::size_t i = 0; i < points.size(); i += 3)
        cgal_points.emplace_back(points[i], points[i + 1], points[i + 2]);

    // new_to_old[new_index] = old_index
    permutation.vertex_new_to_old.resize(points.size() / 3);
    std::iota(permutation.vertex_new_to_old.begin(), permutation.vertex_new_to_old.end(), 0);

    {
        auto point_map = CGAL::make_property_map(cgal_points);

        using Traits =
            CGAL::Spatial_sort_traits_adapter_3<
                Kernel,
                decltype(point_map)>;

        // should move to a concurrency tag template parameter in the future.
        // but currently "incompatible" with openMP in Mesh_smoothing_3
        CGAL::hilbert_sort<CGAL::Sequential_tag>(
            permutation.vertex_new_to_old.begin(),
            permutation.vertex_new_to_old.end(),
            Traits(point_map));
    }

    // Construct inverse permutation: old -> new
    permutation.vertex_old_to_new.resize(points.size() / 3);

    for (std::size_t new_i = 0; new_i < points.size() / 3; ++new_i)
        permutation.vertex_old_to_new[permutation.vertex_new_to_old[new_i]] = new_i;

    // Actually reorder vertices.
    std::vector<double> sorted_points(points.size());

    for (std::size_t new_i = 0; new_i < points.size() / 3; ++new_i) {
        sorted_points[3 * new_i]     = points[3 * permutation.vertex_new_to_old[new_i]];
        sorted_points[3 * new_i + 1] = points[3 * permutation.vertex_new_to_old[new_i] + 1];
        sorted_points[3 * new_i + 2] = points[3 * permutation.vertex_new_to_old[new_i] + 2];
    }

    points.swap(sorted_points);

    // Update tetrahedron connectivity.
    for (auto& tet : tets)
        for (unsigned& v : tet)
            v = permutation.vertex_old_to_new[v];


    // -------------------------------------------------------------------------
    // Tetrahedra
    //
    // Sort according to their barycenters.
    // -------------------------------------------------------------------------

    std::vector<Point_3> tet_centers;
    tet_centers.reserve(tets.size());

    for (const auto& tet : tets)
    {
        double x = (points[3 * tet[0] + 0] + points[3 * tet[1] + 0] + points[3 * tet[2] + 0] + points[3 * tet[3] + 0]) / 4.0;
        double y = (points[3 * tet[0] + 1] + points[3 * tet[1] + 1] + points[3 * tet[2] + 1] + points[3 * tet[3] + 1]) / 4.0;
        double z = (points[3 * tet[0] + 2] + points[3 * tet[1] + 2] + points[3 * tet[2] + 2] + points[3 * tet[3] + 2]) / 4.0;
        tet_centers.emplace_back(x, y, z);
    }

    // new_to_old[new_tet] = old_tet
    permutation.tet_new_to_old.resize(tets.size());
    std::iota(permutation.tet_new_to_old.begin(), permutation.tet_new_to_old.end(), 0);

    {
        auto point_map = CGAL::make_property_map(tet_centers);

        using Traits =
            CGAL::Spatial_sort_traits_adapter_3<
                Kernel,
                decltype(point_map)>;

        CGAL::hilbert_sort<CGAL::Sequential_tag>(
            permutation.tet_new_to_old.begin(),
            permutation.tet_new_to_old.end(),
            Traits(point_map));
    }

    permutation.tet_old_to_new.resize(tets.size());

    for (std::size_t new_i = 0; new_i < tets.size(); ++new_i)
        permutation.tet_old_to_new[permutation.tet_new_to_old[new_i]] = new_i;

    // Actually reorder tetrahedra.
    std::vector<std::array<unsigned, 4>> sorted_tets(tets.size());

    for (std::size_t new_i = 0; new_i < tets.size(); ++new_i)
        sorted_tets[new_i] = tets[permutation.tet_new_to_old[new_i]];

    tets.swap(sorted_tets);

    return permutation;
}

} } // end of CGAL::Mesh_smoothing_3_internal namespace

#endif // CGAL_MESH_SMOOTHING_3_INTERNAL_HILBERT_SORT_H