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

#include <CGAL/Mesh_smoothing_3/internal/type_definitions.h>

#include <CGAL/centroid.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Kernel_traits.h>

#include <array>
#include <iterator>
#include <map>
#include <vector>
#include <tuple>
#include <iostream>

namespace CGAL {

namespace Mesh_smoothing_3 {

/// \cgalAdvancedBegin
/*!
 * \ingroup pkgMeshSmoothing3Projection
 *
 * \brief specifies the weight used for projection onto a tangent space.
 */
enum class Projection_weight_mode
{
    DEFAULT, ///< Standard projection weight (1.).
    STRONG,  ///< Strong projection constraint (10.).
    SOFT,    ///< Weak projection constraint (1e-3).
    NONE,    ///< Disable projection.
    CUSTOM   ///< Use the weight returned by `TangentSpace::custom_weight()`.
};
/// \cgalAdvancedEnd


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
* \brief provides projection to a mesh defined in a tetrahedral mesh model of `MeshComplex_3InTriangulation_3`.
*
* The class creates an `AABB_tree` for each patch and curve on the mesh.
* It then defines queries re-projecting entities on the mesh depending on their patch/curve index.
*
* @tparam C3t3 model of `MeshComplex_3InTriangulation_3`
*
* \cgalModels{ConstructTangentSpace}
*
\sa `CGAL::boundary_aware_mesh_smoothing()`
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
        CConstructor

        \param c3t3 is the mesh used for projection.
    */
    C3t3_mesh_projector(C3t3 const& c3t3)
    : _c3t3(c3t3)
    {
        build_trees();
    }
private:
    C3t3 _c3t3;

    using Facet_tree = typename Mesh_smoothing_3_internal::Facet_tree<C3t3>;
    using Edge_tree = typename Mesh_smoothing_3_internal::Edge_tree<C3t3>;

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
\sa `CGAL::boundary_aware_mesh_smoothing()`
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
        return Tangent_space{Point_3(), typename Tangent_space::Vector_3(), Projection_weight_mode::NONE, 0.};
    }

    Tangent_space curve_edge_projection_line(Curve_edge, std::array<Point_3,2> const &) const {
        return Tangent_space{Point_3(), typename Tangent_space::Vector_3(), Projection_weight_mode::NONE, 0.};
    }

};


/*!
 * \ingroup pkgMeshSmoothing3Projection
 *
 * \brief provides projection functions onto the surface of a polyhedral mesh domain.
 *
 * This class adapts a polyhedral mesh domain into a model of `ConstructTangentSpace`.
 * For each surface facet, the tangent plane is constructed from the closest point
 * and triangle on the polyhedral surface.
 *
 * Curve projection is disabled, as a polyhedral mesh domain without features does
 * not provide geometric information about feature curves.
 *
 * @tparam MeshDomain a model of `MeshDomain_3` providing an AABB tree of its
 * polyhedral boundary, such as `CGAL::Polyhedral_mesh_domain_3`
 *
 * \cgalModels{ConstructTangentSpace}
 *
 * \sa `CGAL::boundary_aware_mesh_smoothing()`
 * \sa `CGAL::Mesh_smoothing_3::Polyhedral_mesh_domain_with_features_projector`
 * \sa `CGAL::Mesh_smoothing_3::C3t3_mesh_projector`
 *
 */
template<typename MeshDomain>
class Polyhedral_mesh_domain_projector
{
public:
    using Point_3 = typename MeshDomain::Point_3;
    using Geom_traits = typename CGAL::Kernel_traits<Point_3>::Kernel;
    using Vector_3 = typename Geom_traits::Vector_3;
    using Tangent_space = Mesh_smoothing_3::Tangent_space<Geom_traits>;

    /*!
     * Constructor
     *
     * \param domain the polyhedral mesh domain used for projection.
     */
    explicit Polyhedral_mesh_domain_projector(const MeshDomain& domain)
        : _domain(domain)
    {}

    template<typename Patch_face>
    Tangent_space
    patch_face_projection_plane(const Patch_face&,
                                const std::vector<Point_3>& face_points) const
    {
        const Point_3 face_center =
            CGAL::centroid(face_points.begin(), face_points.end());

        const auto closest =
            _domain.aabb_tree().closest_point_and_primitive(face_center);

        using AABB_primitive = typename MeshDomain::AABB_primitive;

        const AABB_primitive primitive(
            closest.second.first,
            *closest.second.second);

        const auto triangle = primitive.datum();

        return Tangent_space{
            closest.first,
            CGAL::normal(
                triangle.vertex(0),
                triangle.vertex(1),
                triangle.vertex(2))
        };
    }

    template<typename Curve_edge>
    Tangent_space
    curve_edge_projection_line(const Curve_edge&,
                               const std::array<Point_3, 2>& edge_points) const
    {
        return Tangent_space{
            CGAL::midpoint(edge_points[0], edge_points[1]),
            edge_points[1] - edge_points[0],
            Projection_weight_mode::NONE
        };
    }

protected:
    const MeshDomain& domain() const
    {
        return _domain;
    }

private:
    const MeshDomain& _domain;
};


/*!
 * \ingroup pkgMeshSmoothing3Projection
 *
 * \brief provides projection functions onto the surface and feature curves of
 * a polyhedral mesh domain with features.
 *
 * This class adapts a polyhedral mesh domain with features into a model of
 * `ConstructTangentSpace`. Surface tangent planes are constructed from the
 * polyhedral boundary, as in `Polyhedral_mesh_domain_projector`.
 *
 * For feature edges, the tangent line is constructed from the segment of the
 * corresponding feature polyline that is closest to the edge center.
 *
 * @tparam MeshDomain a model of `MeshDomainWithFeatures_3` providing a
 * polyhedral boundary and feature polylines, such as
 * `CGAL::Polyhedral_mesh_domain_with_features_3`
 *
 * \cgalModels{ConstructTangentSpace}
 *
 * \sa `CGAL::boundary_aware_mesh_smoothing()`
 * \sa `CGAL::Mesh_smoothing_3::Polyhedral_mesh_domain_projector`
 * \sa `CGAL::Mesh_smoothing_3::C3t3_mesh_projector`
 *
 */
template<typename MeshDomain>
class Polyhedral_mesh_domain_with_features_projector
    : public Polyhedral_mesh_domain_projector<MeshDomain>
{
    using Base = Polyhedral_mesh_domain_projector<MeshDomain>;

public:
    using Geom_traits = typename Base::Geom_traits;
    using Point_3 = typename Base::Point_3;
    using Vector_3 = typename Base::Vector_3;
    using Tangent_space = typename Base::Tangent_space;
    using Segment_3 = typename Geom_traits::Segment_3;

    /*!
     * Constructor
     *
     * \param domain the polyhedral mesh domain with features used for projection.
     */
    explicit Polyhedral_mesh_domain_with_features_projector(
        const MeshDomain& domain)
        : Base(domain)
    {}

    template<typename Curve_edge>
    Tangent_space
    curve_edge_projection_line(const Curve_edge& curve_edge,
                               const std::array<Point_3, 2>& edge_points) const
    {
        const Point_3 edge_center =
            CGAL::midpoint(edge_points[0], edge_points[1]);

        // locate_point() returns the source of the polyline segment
        // closest to edge_center on the requested curve.
        const auto segment_source =
            this->domain().locate_point(curve_edge.first, edge_center);

        const Point_3 source = *segment_source;
        const Point_3 target = *std::next(segment_source);

        const Segment_3 segment(source, target);

        const Point_3 projected =
            Geom_traits().construct_projected_point_3_object()(
                segment, edge_center);

        return Tangent_space{
            projected,
            target - source
        };
    }
};

/*!
 * \ingroup pkgMeshSmoothing3Projection
 *
 * \brief provides projection functions onto a surface represented by a
 * signed-distance function.
 *
 * This class adapts a signed-distance function into a model of
 * `ConstructTangentSpace`. The function must return both its signed distance
 * and gradient at a queried point. Negative and positive values represent the
 * two sides of the surface, and the zero level set represents the target
 * surface.
 *
 * Surface points are obtained by iteratively moving along the gradient toward
 * the zero level set. For an exact signed-distance function, whose gradient has
 * unit norm, a single iteration gives the normal projection. Several
 * iterations are supported to accommodate approximate signed-distance fields,
 * for example fields interpolated from a regular grid.
 *
 * Convergence is determined from the magnitude of the projection displacement,
 * rather than from the signed-distance value. This ensures that a large first
 * projection landing exactly on the zero level set is followed by another
 * evaluation before convergence is reported.
 *
 * Curve projection is disabled.
 *
 * @tparam GeomTraits a geometric traits class
 * @tparam Function a callable object taking a `Point_3` and returning
 * `std::pair<FT, Vector_3>`, containing respectively the signed distance
 * and its gradient. This gradient must never be strictly zero.
 *
 * \cgalModels{ConstructTangentSpace}
 *
 * \sa `CGAL::boundary_aware_mesh_smoothing()`
 * \sa `CGAL::Mesh_smoothing_3::Polyhedral_mesh_domain_projector`
 * \sa `CGAL::Mesh_smoothing_3::C3t3_mesh_projector`
 */
template<typename GeomTraits, typename Function>
class Signed_distance_function_projector
{
public:
    using FT = typename GeomTraits::FT;
    using Point_3 = typename GeomTraits::Point_3;
    using Vector_3 = typename GeomTraits::Vector_3;
    using Tangent_space = Mesh_smoothing_3::Tangent_space<GeomTraits>;

    /*!
     * Constructor
     *
     * \param function the signed-distance function and gradient used for
     * projection
     * \param max_projection_iterations maximum number of projection iterations
     * \param tolerance positional tolerance used to stop the projection,
     * expressed in the units of the input geometry
     */
    explicit Signed_distance_function_projector(
        Function function,
        std::size_t max_projection_iterations = 10,
        FT tolerance = FT(1e-8))
        : _function(std::move(function))
        , _max_projection_iterations(max_projection_iterations)
        , _tolerance(tolerance)
    {
        CGAL_precondition(_tolerance >= FT(0));
    }

    template<typename Patch_face>
    Tangent_space
    patch_face_projection_plane(const Patch_face&,
                                const std::vector<Point_3>& face_points) const
    {
        CGAL_precondition(!face_points.empty());

        Point_3 projected =
            CGAL::centroid(face_points.begin(), face_points.end());

        const FT squared_tolerance = _tolerance * _tolerance;

        for(std::size_t i = 0; i < _max_projection_iterations; ++i)
        {
            const auto [distance, gradient] = _function(projected);
            const FT squared_norm = gradient.squared_length();

            CGAL_precondition(squared_norm != FT(0));

            // For an exact signed-distance function, ||gradient|| = 1 and
            // this reduces to: displacement = distance * gradient.
            //
            // Keeping the normalization makes the projection more robust to
            // approximate signed-distance fields.
            const Vector_3 displacement =
                (distance / squared_norm) * gradient;

            projected = projected - displacement;

            // Use positional convergence instead of |distance|. In
            // particular, a large first displacement that happens to land
            // exactly on the zero level set does not immediately terminate
            // the iteration.
            if(displacement.squared_length() <= squared_tolerance)
                break;
        }

        // Re-evaluate at the projected position so that the returned tangent
        // plane corresponds to the final surface point.
        const auto [distance, gradient] = _function(projected);
        CGAL_USE(distance);

        CGAL_precondition(gradient.squared_length() != FT(0));

        return Tangent_space{projected, gradient};
    }

    template<typename Curve_edge>
    Tangent_space
    curve_edge_projection_line(const Curve_edge&,
                               const std::array<Point_3, 2>& edge_points) const
    {
        return Tangent_space{
            CGAL::midpoint(edge_points[0], edge_points[1]),
            edge_points[1] - edge_points[0],
            Projection_weight_mode::NONE
        };
    }

private:
    std::decay_t<Function> _function;
    std::size_t _max_projection_iterations;
    FT _tolerance;
};

} // namespace Mesh_smoothing_3

}

#endif // CGAL_MESH_SMOOTHING_3_PROJECTORS_H
