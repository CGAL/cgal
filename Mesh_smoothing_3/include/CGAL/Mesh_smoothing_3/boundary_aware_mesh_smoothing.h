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

#ifndef CGAL_MESH_SMOOTHING_3_BOUNDARY_AWARE_MESH_SMOOTHING_H
#define CGAL_MESH_SMOOTHING_3_BOUNDARY_AWARE_MESH_SMOOTHING_H

#include <CGAL/license/Mesh_smoothing_3.h>

#include <CGAL/Mesh_smoothing_3/Mesh_smoothing_3.h>

#include <CGAL/Named_function_parameters.h>


namespace CGAL {

/*!
* \ingroup pkgMeshSmoothing3Fonctions
* Smooth a tetrahedral mesh while preserving the boundary and curve features.
*
* This function takes as input a `Mesh_complex_3_in_triangulation_3` and will iteratively
* move its vertices to improve the quality of the tetrahedra while preserving (softly) the boundary and curve features.
* Vertices will stay close to their original entity (surface patch or curve). Corners are locked on their original position.
* The code is capable of recovering curvature discontinuities (sharp features) in a patch, 
* and is capable of recovering inverted tetrahedra in the mesh. 
* A valid input mesh is guaranteed to remain valid after and during the smoothing.
*
* @tparam C3t3 is a `Mesh_complex_3_in_triangulation_3` type.
* @tparam NamedParameters is a sequence of named parameters.    
* 
* @param c3t3 is the tetrahedral mesh to smooth.
* @param np is an optional sequence of \ref bgl_namedparameters "Named Parameters"
*          among the ones listed below
* \cgalNamedParamsBegin
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{the number of iterations for the full sequence of atomic operations
*                           performed (listed in the above description)}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`100`}
*   \cgalParamNEnd
*   \cgalParamNBegin{verbose}
*     \cgalParamDescription{If `true`, the function prints information about the smoothing
*                           process to the standard output.}
*     \cgalParamType{`bool`}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template<typename C3t3,
         typename NamedParameters = parameters::Default_named_parameters>
void boundary_aware_mesh_smoothing  (
    C3t3& c3t3,
    NamedParameters const & np = parameters::default_values()
) 
{

    using parameters::choose_parameter;
    using parameters::get_parameter;

    bool verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);
    std::size_t max_iterations = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 100);

    C3t3 ref (static_cast<const C3t3>(c3t3));
    
    // AABB trees of the surface and curve features of the c3t3
    std::map<C3t3::Surface_patch_index, std::vector<C3t3::Facet>> facets_by_patch;
    for (auto f : ref.facets_in_complex())
    {
        facets_by_patch[ref.surface_patch_index(f)].push_back(f);
    }
    
    std::map<C3t3::Curve_index, std::vector<C3t3::Edge>> edges_by_curve;
    for (auto e : ref.edges_in_complex())
    {
        edges_by_curve[ref.curve_index(e)].push_back(e);
    }

    // AABB traits for C3t3 facets and edges
    using K = typename C3t3::Triangulation::Geom_traits;
    using Tr = typename C3t3::Triangulation;
    struct Facet_to_triangle_property_map
    {
        typedef C3t3::Facet key_type;
        typedef K::Triangle_3 value_type;
        typedef value_type reference;
        typedef boost::readable_property_map_tag category;
    };

    inline Facet_to_triangle_property_map::reference
    get(const Facet_to_triangle_property_map&, const C3t3::Facet& f)
    {
        const Tr::Cell_handle cell = f.first;
        const int i = f.second;
        const K::Point_3& p0 = cell->vertex((i + 1) % 4)->point().point();
        const K::Point_3& p1 = cell->vertex((i + 2) % 4)->point().point();
        const K::Point_3& p2 = cell->vertex((i + 3) % 4)->point().point();
        return K::Triangle_3(p0, p1, p2);
    }

    struct Facet_to_point_property_map
    {
        typedef C3t3::Facet key_type;
        typedef K::Point_3 value_type;
        typedef value_type reference;
        typedef boost::readable_property_map_tag category;
    };

    inline Facet_to_point_property_map::reference
    get(const Facet_to_point_property_map&, const C3t3::Facet& f)
    {
        const Tr::Cell_handle cell = f.first;
        const int i = f.second;
        return cell->vertex((i + 1) % 4)->point().point();
    }

    struct Edge_to_segment_property_map
    {
        typedef C3t3::Edge key_type;
        typedef K::Segment_3 value_type;
        typedef value_type reference;
        typedef boost::readable_property_map_tag category;
    };

    inline Edge_to_segment_property_map::reference
    get(const Edge_to_segment_property_map&, const C3t3::Edge& e)
    {
        const Tr::Cell_handle cell = e.first;
        return K::Segment_3(cell->vertex(e.second)->point().point(), cell->vertex(e.third)->point().point());
    }

    struct Edge_to_point_property_map
    {
        typedef C3t3::Edge key_type;
        typedef K::Point_3 value_type;
        typedef value_type reference;
        typedef boost::readable_property_map_tag category;
    };

    inline Edge_to_point_property_map::reference
    get(const Edge_to_point_property_map&, const C3t3::Edge& e)
    {
        const Tr::Cell_handle cell = e.first;
        return cell->vertex(e.second)->point().point();
    }

    typedef CGAL::AABB_primitive<C3t3::Facet,
                                Facet_to_triangle_property_map,
                                Facet_to_point_property_map,
                                CGAL::Tag_false,
                                CGAL::Tag_false> Facet_primitive;
    typedef CGAL::AABB_traits_3<K, Facet_primitive> Facet_traits;
    typedef CGAL::AABB_tree<Facet_traits> Facet_tree;

    typedef CGAL::AABB_primitive<C3t3::Edge,
                                Edge_to_segment_property_map,
                                Edge_to_point_property_map,
                                CGAL::Tag_false,
                                CGAL::Tag_false> Edge_primitive;
    typedef CGAL::AABB_traits_3<K, Edge_primitive> Edge_traits;
    typedef CGAL::AABB_tree<Edge_traits> Edge_tree;

    std::map<C3t3::Surface_patch_index, Facet_tree> facet_trees;
    for (auto& kv : facets_by_patch)
    {
        facet_trees.try_emplace(kv.first, kv.second.begin(), kv.second.end());
    }

    std::map<C3t3::Curve_index, Edge_tree> edge_trees;
    for (auto& kv : edges_by_curve)
    {
        edge_trees.try_emplace(kv.first, kv.second.begin(), kv.second.end());
    }

    if (verbose) std::cout << "Built " << facet_trees.size() << " facet AABB trees and " << edge_trees.size() << " edge AABB trees." << std::endl;


    // Smoothing

    CGAL::Mesh_smoothing_3::C3t3_smoother smoother(c3t3);

    for (auto c : c3t3.vertices_in_complex())
    {
        smoother.set_vertex_lock(c, true);
    }

    smoother.set_boundary_query([&](K::Point_3 const &pt, C3t3::Surface_patch_index id, double radius) -> std::tuple<K::Point_3, K::Vector_3, double> {
        auto res = facet_trees.at(id).closest_point_and_primitive(pt);
        K::Point_3 closest_point = res.first;
        const auto triangle = ref.triangulation().triangle(res.second);
        K::Vector_3 normal = CGAL::unit_normal(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));
        return std::make_tuple(closest_point, normal, 1.);
    });

    smoother.set_curves_query([&](K::Point_3 const &pt, C3t3::Curve_index id, double radius) -> std::tuple<K::Point_3, K::Vector_3, double> {
        auto res = edge_trees.at(id).closest_point_and_primitive(pt);
        K::Point_3 closest_point = res.first;
        const auto segment = ref.triangulation().segment(res.second);
        K::Vector_3 direction = (segment.target() - segment.source());
        if (direction.squared_length() > 1e-8)
        {
            direction /= CGAL::sqrt(direction.squared_length());
        }
        return std::make_tuple(closest_point, direction, 1.);
    });

    smoother.set_verbose(verbose);
    smoother.set_max_number_of_iteration(max_iterations);
    smoother.run();
}

} 

#endif // CGAL_MESH_SMOOTHING_3_BOUNDARY_AWARE_MESH_SMOOTHING_H

