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

public:
    Tangent_space<Geom_traits> patch_face_projection_plane(Patch_face patch_face, std::vector<Point_3> const &face_points) const {
        Tangent_space<Geom_traits> projection;
        Point_3 face_center = CGAL::centroid(face_points.begin(), face_points.end());
        auto res = _facet_trees.at(patch_face.first).closest_point_and_primitive(face_center);
        projection._origin = res.first;
        const auto triangle = _c3t3.triangulation().triangle(res.second);
        projection._vector = CGAL::unit_normal(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));
        return projection;
    }

    Tangent_space<Geom_traits> curve_edge_projection_line(Curve_edge curve_edge, std::array<Point_3,2> const &edge_points) const {
        Tangent_space<Geom_traits> projection;
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

public:

    Tangent_space<Geom_traits> patch_face_projection_plane(Patch_face, std::vector<Point_3> const &) const {
        return Tangent_space<Geom_traits>{Point_3(), Geom_traits::Vector_3(), Projection_weight_mode::NONE, 0.};
    }

    Tangent_space<Geom_traits> curve_edge_projection_line(Curve_edge, std::array<Point_3,2> const &) const {
        return Tangent_space<Geom_traits>{Point_3(), Geom_traits::Vector_3(), Projection_weight_mode::NONE, 0.};
    }

};

} // namespace Mesh_smoothing_3

/*!
* \ingroup pkgMeshSmoothing3Functions
* Smooth a tetrahedral mesh while preserving the boundary and curve features.
*
* This function takes as input a `Mesh_complex_3_in_triangulation_3` and will iteratively
* move its vertices to improve the quality of the tetrahedra while preserving (softly) the boundary and curve features through projection queries.
* Vertices will stay close to their original entity (surface patch or curve). Corners are locked on their original position.
* The algorithm will recover curvature discontinuities (sharp features) in a patch,
* and is capable of recovering inverted tetrahedra in the mesh.
* An input mesh with only positively oriented elements is guaranteed to remain so during the smoothing procedure.
*
* \warning This function updates vertex positions without modifying connectivity, and therefore does not preserve Delaunay properties.
*
* @tparam C3t3 model of `MeshComplex_3InTriangulation_3`
* @tparam CTS model of `ConstructTangentSpace`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
*
* @param c3t3 tetrahedral mesh to smooth.
* @param cts support to re-project faces and edges on their respective patches and curves.
* @param np optional sequence of \ref bgl_namedparameters "Named Parameters"
*          among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{Maximum nb of iterations of the smoothing algorithm .
                            Algorithm will stop before if it reaches convergence.
                            Untangling usually requires more iterations (up to thousands) for hard cases. }
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`100`}
*   \cgalParamNEnd
*   \cgalParamNBegin{verbose}
*     \cgalParamDescription{If `true`, the function prints information about the smoothing
*                           process to the standard output.}
*     \cgalParamType{`bool`}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
*   \cgalParamNBegin{vertex_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each vertex of the triangulation in the `c3t3`.}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `C3t3::Vertex_handle`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no vertex is constrained}
*     \cgalParamExtra{A constrained vertex will not be moved by smoothing.}
*   \cgalParamNEnd
*   \cgalParamNBegin{edge_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of the triangulation in the `c3t3`.}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `std::pair<C3t3::Vertex_handle, C3t3::Vertex_handle>`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no edge is constrained}
*     \cgalParamExtra{Vertices of a constrained edge will not be moved by the smoothing.}
*   \cgalParamNEnd
*   \cgalParamNBegin{facet_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each facet of the triangulation in the `c3t3`.}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `C3t3::Facet`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no facet is constrained}
*     \cgalParamExtra{Vertices of a constrained facet will not be moved by smoothing.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template<typename C3t3,
         typename CTS,
         typename NamedParameters = parameters::Default_named_parameters>
void boundary_aware_mesh_smoothing  (
    C3t3& c3t3,
    CTS const & cts,
    NamedParameters const & np = parameters::default_values()
)
{
    using parameters::choose_parameter;
    using parameters::get_parameter;

    const bool verbose =
        choose_parameter(get_parameter(np, internal_np::verbose), false);
    const std::size_t max_iterations =
        choose_parameter(get_parameter(np, internal_np::number_of_iterations), 100u);


    typedef typename internal_np::Lookup_named_param_def <
        internal_np::vertex_is_constrained_t,
        NamedParameters,
        Static_boolean_property_map<typename C3t3::Vertex_handle, false> // default (no constraint pmap)
    > ::type VCMap;
    VCMap vcmap = choose_parameter<Static_boolean_property_map<typename C3t3::Vertex_handle, false>>(get_parameter(np, internal_np::vertex_is_constrained));

    typedef typename internal_np::Lookup_named_param_def <
        internal_np::edge_is_constrained_t,
        NamedParameters,
        Static_boolean_property_map<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>, false> // default (no constraint pmap)
    > ::type ECMap;
    ECMap ecmap = choose_parameter<Static_boolean_property_map<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>, false>>(get_parameter(np, internal_np::edge_is_constrained));

    typedef typename internal_np::Lookup_named_param_def <
        internal_np::facet_is_constrained_t,
        NamedParameters,
        Static_boolean_property_map<typename C3t3::Facet, false> // default (no constraint pmap)
    > ::type FCMap;
    FCMap fcmap = choose_parameter<Static_boolean_property_map<typename C3t3::Facet, false>>(get_parameter(np, internal_np::facet_is_constrained));


    // C3t3_smoother will automatically mark facets and features
    CGAL::Mesh_smoothing_3::C3t3_smoother smoother(c3t3);

    // locks through property maps
    for (auto v : c3t3.triangulation().finite_vertex_handles()) {
        if (!get(vcmap, v)) continue;
        smoother.set_vertex_lock(v, true);
    }

    for (const auto& e : c3t3.triangulation().finite_edges()) {
        const auto vertices = c3t3.triangulation().vertices(e);
        if (get(ecmap, {vertices[0], vertices[1]}) || get(ecmap, {vertices[1], vertices[0]})) {
            smoother.set_vertex_lock(vertices[0], true);
            smoother.set_vertex_lock(vertices[1], true);
        }
    }

    for (const auto& f : c3t3.triangulation().finite_facets()) {
        if (get(fcmap, f)) {
            for (const auto& v : c3t3.triangulation().vertices(f))
                smoother.set_vertex_lock(v, true);
        }
    }

    // locks corners
    for (auto c : c3t3.vertices_in_complex()) {
        smoother.set_vertex_lock(c, true);
    }

    using Point_3 = typename  C3t3::Triangulation::Geom_traits::Point_3;
    using Tangent_space = Mesh_smoothing_3::Tangent_space<typename C3t3::Triangulation::Geom_traits>;

    auto proj_weight = [](Tangent_space const &ts) {
        using Mesh_smoothing_3::Projection_weight_mode;
        switch (ts.projection_mode())
        {
        case Projection_weight_mode::DEFAULT:
            return 1.;
            break;
        case Projection_weight_mode::STRONG:
            return 10.;
            break;
        case Projection_weight_mode::SOFT:
            return 1e-3;
            break;
        case Projection_weight_mode::NONE:
            return 0.;
            break;
        case Projection_weight_mode::CUSTOM:
            return ts.custom_weight();
            break;
        default:
            return 1.;
            break;
        }
    };

    // converting Projector to Mesh_smoothing_3 queries
    smoother.set_boundary_query([&](std::vector<Point_3> const& pts, std::pair<typename C3t3::Surface_patch_index, typename C3t3::Facet> patch_face) {
        Tangent_space proj = cts.patch_face_projection_plane(patch_face, pts);
        return std::make_tuple(proj.origin(), proj.vector(), proj_weight(proj));
    });

    smoother.set_curves_query([&](std::array<Point_3, 2> const& pts, std::pair<typename C3t3::Curve_index, typename C3t3::Edge> curve_edge) {
        Tangent_space proj = cts.curve_edge_projection_line(curve_edge, pts);
        return std::make_tuple(proj.origin(), proj.vector(), proj_weight(proj));
    });

    smoother.set_verbose(verbose);
    smoother.set_max_number_of_iteration(max_iterations);
    smoother.run();
}

}

#endif // CGAL_MESH_SMOOTHING_3_BOUNDARY_AWARE_MESH_SMOOTHING_H

