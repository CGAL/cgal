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

#ifndef CGAL_MESH_SMOOTHING_3_PROJECTORS_H
#define CGAL_MESH_SMOOTHING_3_PROJECTORS_H


#include <CGAL/license/Mesh_smoothing_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_primitive.h>
#include <CGAL/centroid.h>

#include <map>
#include <vector>
#include <tuple>
#include <iostream>

namespace CGAL {

namespace Mesh_smoothing_3 {
namespace internal {

using CGAL::Mesh_smoothing_3::cgal_types::get_point;

template<typename C3t3>
struct Facet_to_triangle_property_map
{
    using key_type = typename C3t3::Facet;
    using K = typename C3t3::Triangulation::Geom_traits;
    using value_type = typename K::Triangle_3;
    using reference = value_type;
    using category = boost::readable_property_map_tag;
};

template<typename C3t3>
inline typename Facet_to_triangle_property_map<C3t3>::reference
get(const Facet_to_triangle_property_map<C3t3>&, const typename C3t3::Facet& f)
{
    using Tr = typename C3t3::Triangulation;
    using K = typename Tr::Geom_traits;
    const typename Tr::Cell_handle cell = f.first;
    const int i = f.second;
    const typename K::Point_3& p0 = get_point<C3t3>(cell->vertex((i + 1) % 4));
    const typename K::Point_3& p1 = get_point<C3t3>(cell->vertex((i + 2) % 4));
    const typename K::Point_3& p2 = get_point<C3t3>(cell->vertex((i + 3) % 4));
    return typename K::Triangle_3(p0, p1, p2);
}

template<typename C3t3>
struct Facet_to_point_property_map
{
    using key_type = typename C3t3::Facet;
    using K = typename C3t3::Triangulation::Geom_traits;
    using value_type = typename K::Point_3;
    using reference = value_type;
    using category = boost::readable_property_map_tag;
};

template<typename C3t3>
inline typename Facet_to_point_property_map<C3t3>::reference
get(const Facet_to_point_property_map<C3t3>&, const typename C3t3::Facet& f)
{
    using Tr = typename C3t3::Triangulation;
    const typename Tr::Cell_handle cell = f.first;
    const int i = f.second;
    return get_point<C3t3>(cell->vertex((i + 1) % 4));
}

template<typename C3t3>
struct Edge_to_segment_property_map
{
    using key_type = typename C3t3::Edge;
    using K = typename C3t3::Triangulation::Geom_traits;
    using value_type = typename K::Segment_3;
    using reference = value_type;
    using category = boost::readable_property_map_tag;
};

template<typename C3t3>
inline typename Edge_to_segment_property_map<C3t3>::reference
get(const Edge_to_segment_property_map<C3t3>&, const typename C3t3::Edge& e)
{
    using Tr = typename C3t3::Triangulation;
    using K = typename Tr::Geom_traits;
    const typename Tr::Cell_handle cell = e.first;
    return typename K::Segment_3(get_point<C3t3>(cell->vertex(e.second)),
                                 get_point<C3t3>(cell->vertex(e.third)));
}

template<typename C3t3>
struct Edge_to_point_property_map
{
    using key_type = typename C3t3::Edge;
    using K = typename C3t3::Triangulation::Geom_traits;
    using value_type = typename K::Point_3;
    using reference = value_type;
    using category = boost::readable_property_map_tag;
};

template<typename C3t3>
inline typename Edge_to_point_property_map<C3t3>::reference
get(const Edge_to_point_property_map<C3t3>&, const typename C3t3::Edge& e)
{
    using Tr = typename C3t3::Triangulation;
    const typename Tr::Cell_handle cell = e.first;
    return get_point<C3t3>(cell->vertex(e.second));
}

} // namespace internal


/*!
 * \ingroup pkgMeshSmoothing3Projection
 *
 * \brief Specifies the weight used for projection onto a tangent space.
 */
enum class Projection_weight_mode
{
    DEFAULT, ///< Standard projection weight (1.).
    STRONG,  ///< Strong projection constraint (10.).
    SOFT,    ///< Weak projection constraint (1e-3).
    NONE,    ///< Disable projection.
    CUSTOM   ///< Use the weight returned by `TangentSpace::custom_weight()`.
};

/* not documented but:
* \cgalModels{TangentSpace}
*/
template <typename GeomTraits>
struct Tangent_space {
    using Point_3 = typename GeomTraits::Point_3;
    using Vector_3 = typename GeomTraits::Vector_3;

    Point_3 _origin = Point_3();
    Vector_3 _vector = Vector_3();
    Projection_weight_mode _mode = Projection_weight_mode::DEFAULT;
    double _weight = 1.;

    Point_3 origin() const { return _origin; }
    Vector_3 vector() const { return _vector; }
    Projection_weight_mode projection_mode() const { return _mode; }
    double custom_weight() const { return _weight; }

};

/*!
* \ingroup pkgMeshSmoothing3Projection
*
* \brief provides projection functions to a mesh defined in a `Mesh_complex_3_in_triangulation_3`
*
* The class `Mesh_smoother` creates an AABB_tree for each patch and curve on the mesh.
* It then defines queries re-projecting entities on the mesh depending on their patch/curve index.
*
* @tparam C3t3 model of `MeshComplex_3InTriangulation_3`
*
* \cgalModels{ConstructTangentSpace}
*
\sa `CGAL::boundary_aware_mesh_smoothing`
\sa `CGAL::Mesh_smoothing_3::C3t3_no_projection`
*
*/
template<typename C3t3>
class C3t3_mesh_projector {
public:
    using Geom_traits = typename C3t3::Triangulation::Geom_traits;
    using Point_3 = typename Geom_traits::Point_3;

    using Facet = typename C3t3::Facet;
    using Surface_patch_index = typename C3t3::Surface_patch_index;
    using Edge = typename C3t3::Edge;
    using Curve_index = typename C3t3::Curve_index;

    using Patch_face = std::pair<Surface_patch_index, Facet>;
    using Curve_edge = std::pair<Curve_index, Edge>;

    using Tangent_space = typename Mesh_smoothing_3::Tangent_space<Geom_traits>;

public:
    Tangent_space patch_face_projection_plane(Patch_face patch_face, std::vector<Point_3> const &face_points) const {
        Tangent_space projection;
        Point_3 face_center = CGAL::centroid(face_points.begin(), face_points.end());
        auto res = _facet_trees.at(patch_face.first).closest_point_and_primitive(face_center);
        projection._origin = res.first;
        const auto triangle = _c3t3.triangulation().triangle(res.second);
        projection._vector = CGAL::unit_normal(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));
        return projection;
    }

    Tangent_space curve_edge_projection_line(Curve_edge curve_edge, std::array<Point_3,2> const &edge_points) const {
        Tangent_space projection;
        Point_3 edge_center = CGAL::midpoint(edge_points[0], edge_points[1]);
        auto res = _edge_trees.at(curve_edge.first).closest_point_and_primitive(edge_center);
        projection._origin =  res.first;
        const auto segment = _c3t3.triangulation().segment(res.second);
        projection._vector = (segment.target() - segment.source());
        if (projection._vector.squared_length() > 1e-8) {
            projection._vector /= CGAL::sqrt(projection._vector.squared_length());
        }
        return projection;
    }


public:

    /*!
        Class constructor

        \param c3t3 contains the mesh used for projection. A copy is done at construction.

    */
    C3t3_mesh_projector(C3t3 const& c3t3)
    : _c3t3(c3t3)
    {
        build_trees();
    }
private:
    C3t3 _c3t3;


    using K = typename C3t3::Triangulation::Geom_traits;

    using Facet_to_triangle_map = internal::Facet_to_triangle_property_map<C3t3>;
    using Facet_to_point_map = internal::Facet_to_point_property_map<C3t3>;
    using Edge_to_segment_map = internal::Edge_to_segment_property_map<C3t3>;
    using Edge_to_point_map = internal::Edge_to_point_property_map<C3t3>;

    using Facet_primitive = CGAL::AABB_primitive<Facet,
                                                Facet_to_triangle_map,
                                                Facet_to_point_map,
                                                CGAL::Tag_false,
                                                CGAL::Tag_false>;
    using Facet_traits = CGAL::AABB_traits_3<K, Facet_primitive>;
    using Facet_tree = CGAL::AABB_tree<Facet_traits>;

    using Edge_primitive = CGAL::AABB_primitive<Edge,
                                                Edge_to_segment_map,
                                                Edge_to_point_map,
                                                CGAL::Tag_false,
                                                CGAL::Tag_false>;
    using Edge_traits = CGAL::AABB_traits_3<K, Edge_primitive>;
    using Edge_tree = CGAL::AABB_tree<Edge_traits>;

    std::map<Surface_patch_index, Facet_tree> _facet_trees;
    std::map<Curve_index, Edge_tree> _edge_trees;

    void build_trees() {
        std::map<Surface_patch_index, std::vector<Facet>> facets_by_patch;
        for (const auto& f : _c3t3.facets_in_complex()) {
            facets_by_patch[_c3t3.surface_patch_index(f)].push_back(f);
        }

        std::map<Curve_index, std::vector<Edge>> edges_by_curve;
        for (const auto& e : _c3t3.edges_in_complex()) {
            edges_by_curve[_c3t3.curve_index(e)].push_back(e);
        }

        for (auto& kv : facets_by_patch) {
            _facet_trees.try_emplace(kv.first, kv.second.begin(), kv.second.end());
        }

        for (auto& kv : edges_by_curve) {
            _edge_trees.try_emplace(kv.first, kv.second.begin(), kv.second.end());
        }

    }
};

/*!
* \ingroup pkgMeshSmoothing3Projection
*
* \brief provides an empty class to disable projections, meaning free boundaries.
*
* @tparam C3t3 model of `MeshComplex_3InTriangulation_3`
*
* \cgalModels{ConstructTangentSpace}
*
\sa `CGAL::boundary_aware_mesh_smoothing`
\sa `CGAL::Mesh_smoothing_3::C3t3_mesh_projector`
*
*/
template<typename C3t3>
class C3t3_no_projection {
public:
    using Geom_traits = typename C3t3::Triangulation::Geom_traits;
    using Point_3 = typename Geom_traits::Point_3;
    using Patch_face = std::pair<typename C3t3::Surface_patch_index, typename C3t3::Facet>;
    using Curve_edge = std::pair<typename C3t3::Curve_index, typename C3t3::Edge>;

    using Tangent_space = Mesh_smoothing_3::Tangent_space<Geom_traits>;

public:

    Tangent_space patch_face_projection_plane(Patch_face, std::vector<Point_3> const &) const {
        return Tangent_space{Point_3(), Geom_traits::Vector_3(), Projection_weight_mode::NONE, 0.};
    }

    Tangent_space curve_edge_projection_line(Curve_edge, std::array<Point_3,2> const &) const {
        return Tangent_space{Point_3(), Geom_traits::Vector_3(), Projection_weight_mode::NONE, 0.};
    }

};

} // namespace Mesh_smoothing_3

}

#endif // CGAL_MESH_SMOOTHING_3_PROJECTORS_H

