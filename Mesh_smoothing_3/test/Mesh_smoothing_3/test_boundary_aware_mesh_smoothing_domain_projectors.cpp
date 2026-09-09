// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/remove_far_points_in_mesh_3.h>

#include <CGAL/Mesh_smoothing_3/boundary_aware_mesh_smoothing.h>
#include <CGAL/Mesh_smoothing_3/projectors.h>

#include <CGAL/centroid.h>
#include <CGAL/config.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <limits>
#include <utility>
#include <vector>

namespace {

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using GT = CGAL::Mesh_3::Robust_intersection_traits_3<K>;

using FT = K::FT;
using Point_3 = K::Point_3;
using Vector_3 = K::Vector_3;
using Segment_3 = K::Segment_3;

using Projection_mode =
    CGAL::Mesh_smoothing_3::Projection_weight_mode;

constexpr double SURFACE_TOLERANCE = 0.25;
constexpr double CURVE_TOLERANCE = 0.25;
constexpr double IMPLICIT_TOLERANCE = 0.20;


// -----------------------------------------------------------------------------
// Common helpers
// -----------------------------------------------------------------------------

bool finite_point(const Point_3& p)
{
    return std::isfinite(CGAL::to_double(p.x())) &&
           std::isfinite(CGAL::to_double(p.y())) &&
           std::isfinite(CGAL::to_double(p.z()));
}

template<typename C3t3>
void assert_all_finite(const C3t3& c3t3)
{
    const auto& tr = c3t3.triangulation();
    const auto construct_point = tr.geom_traits().construct_point_3_object();

    for(const auto v : tr.finite_vertex_handles())
    {
        const Point_3 p = construct_point(v->point());
        assert(finite_point(p));
    }
}

template<typename C3t3, typename MeshDomain>
double max_surface_distance(
    const C3t3& c3t3,
    const MeshDomain& domain)
{
    const auto& tr = c3t3.triangulation();

    double max_distance = 0.;

    for(const auto& facet : c3t3.facets_in_complex())
    {
        const Point_3 center = CGAL::centroid(tr.triangle(facet));

        const double distance = std::sqrt(
            CGAL::to_double(
                domain.aabb_tree().squared_distance(center)));

        max_distance = (std::max)(max_distance, distance);
    }

    return max_distance;
}

const std::array<Segment_3, 12>& cube_edges()
{
    static const std::array<Segment_3, 12> edges = {
        Segment_3(Point_3(-1, -1, -1), Point_3( 1, -1, -1)),
        Segment_3(Point_3(-1,  1, -1), Point_3( 1,  1, -1)),
        Segment_3(Point_3(-1, -1,  1), Point_3( 1, -1,  1)),
        Segment_3(Point_3(-1,  1,  1), Point_3( 1,  1,  1)),

        Segment_3(Point_3(-1, -1, -1), Point_3(-1,  1, -1)),
        Segment_3(Point_3( 1, -1, -1), Point_3( 1,  1, -1)),
        Segment_3(Point_3(-1, -1,  1), Point_3(-1,  1,  1)),
        Segment_3(Point_3( 1, -1,  1), Point_3( 1,  1,  1)),

        Segment_3(Point_3(-1, -1, -1), Point_3(-1, -1,  1)),
        Segment_3(Point_3( 1, -1, -1), Point_3( 1, -1,  1)),
        Segment_3(Point_3(-1,  1, -1), Point_3(-1,  1,  1)),
        Segment_3(Point_3( 1,  1, -1), Point_3( 1,  1,  1))
    };

    return edges;
}

double distance_to_cube_edges(const Point_3& p)
{
    double min_squared_distance =
        (std::numeric_limits<double>::max)();

    for(const Segment_3& edge : cube_edges())
    {
        min_squared_distance = (std::min)(
            min_squared_distance,
            CGAL::to_double(CGAL::squared_distance(p, edge)));
    }

    return std::sqrt(min_squared_distance);
}

template<typename C3t3>
double max_curve_distance_to_cube(const C3t3& c3t3)
{
    const auto& tr = c3t3.triangulation();

    double max_distance = 0.;
    std::size_t nb_curves = 0;

    for(const auto& edge : c3t3.edges_in_complex())
    {
        const Segment_3 segment = tr.segment(edge);
        const Point_3 center =
            CGAL::midpoint(segment.source(), segment.target());

        max_distance = (std::max)(
            max_distance,
            distance_to_cube_edges(center));

        ++nb_curves;
    }

    assert(nb_curves > 0);
    return max_distance;
}

using Polyhedron = typename CGAL::Mesh_polyhedron_3<GT>::type;

Polyhedron load_cube()
{
    std::ifstream input(CGAL::data_file_path("meshes/cube.off"));
    assert(input);

    Polyhedron polyhedron;
    input >> polyhedron;

    assert(input);
    return polyhedron;
}


// -----------------------------------------------------------------------------
// Polyhedral_mesh_domain_projector
// -----------------------------------------------------------------------------

void test_polyhedral_mesh_domain_projector()
{
    using Domain =
        CGAL::Polyhedral_mesh_domain_3<Polyhedron, GT>;

    using Tr =
        typename CGAL::Mesh_triangulation_3<
            Domain,
            typename CGAL::Kernel_traits<Domain>::Kernel>::type;

    using C3t3 =
        CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

    using Criteria = CGAL::Mesh_criteria_3<Tr>;
    using Facet_criteria = typename Criteria::Facet_criteria;
    using Cell_criteria = typename Criteria::Cell_criteria;

    using Projector =
        CGAL::Mesh_smoothing_3::Polyhedral_mesh_domain_projector<
            Domain>;

    Polyhedron polyhedron = load_cube();
    Domain domain(polyhedron);

    Facet_criteria facet_criteria(30., 0.3, 0.03);
    Cell_criteria cell_criteria(3., 0.4);
    Criteria criteria(facet_criteria, cell_criteria);

    C3t3 c3t3 =
        CGAL::make_mesh_3<C3t3>(
            domain,
            criteria,
            CGAL::parameters::no_exude().no_perturb());

    CGAL::remove_far_points_in_mesh_3(c3t3);

    assert(c3t3.number_of_facets_in_complex() > 0);

    Projector projector(domain);

    // Basic direct surface query:
    // only check that the result is finite, non-degenerate, and close
    // to the target surface.
    {
        const auto facet = *c3t3.facets_in_complex().begin();
        const auto triangle = c3t3.triangulation().triangle(facet);

        std::vector<Point_3> points = {
            triangle[0],
            triangle[1],
            triangle[2]
        };

        const auto tangent_space =
            projector.patch_face_projection_plane(
                std::make_pair(
                    c3t3.surface_patch_index(facet),
                    facet),
                points);

        assert(finite_point(tangent_space.origin()));
        assert(tangent_space.vector().squared_length() > FT(0));

        const double distance = std::sqrt(
            CGAL::to_double(
                domain.aabb_tree().squared_distance(
                    tangent_space.origin())));

        assert(distance < SURFACE_TOLERANCE);
    }

    // A domain without features must explicitly disable curve projection.
    {
        const std::array<Point_3, 2> points = {
            Point_3(-1., -1., 0.),
            Point_3(-1.,  1., 0.)
        };

        const auto tangent_space =
            projector.curve_edge_projection_line(0, points);

        assert(
            tangent_space.projection_mode() ==
            Projection_mode::NONE);

        assert(tangent_space.vector().squared_length() > FT(0));
    }

    const auto status =
        CGAL::boundary_aware_mesh_smoothing(
            c3t3,
            projector);

    assert_all_finite(c3t3);
    assert(status.valid_mesh());
    assert(status.nb_invalid_elements == 0);

    // The smoother balances fitting and element quality, so do not require
    // exact reprojection onto the cube.
    assert(
        max_surface_distance(c3t3, domain) <
        SURFACE_TOLERANCE);
}


// -----------------------------------------------------------------------------
// Polyhedral_mesh_domain_with_features_projector
// -----------------------------------------------------------------------------

void test_polyhedral_mesh_domain_with_features_projector()
{
    using Domain =
        CGAL::Polyhedral_mesh_domain_with_features_3<
            GT,
            Polyhedron>;

    using Tr =
        typename CGAL::Mesh_triangulation_3<
            Domain,
            typename CGAL::Kernel_traits<Domain>::Kernel>::type;

    using C3t3 =
        CGAL::Mesh_complex_3_in_triangulation_3<
            Tr,
            typename Domain::Corner_index,
            typename Domain::Curve_index>;

    using Criteria = CGAL::Mesh_criteria_3<Tr>;
    using Edge_criteria = typename Criteria::Edge_criteria;
    using Facet_criteria = typename Criteria::Facet_criteria;
    using Cell_criteria = typename Criteria::Cell_criteria;

    using Projector =
        CGAL::Mesh_smoothing_3::
            Polyhedral_mesh_domain_with_features_projector<Domain>;

    Polyhedron polyhedron = load_cube();

    Domain domain(
        polyhedron,
        &CGAL::get_default_random());

    domain.detect_features();

    Edge_criteria edge_criteria(0.3);
    Facet_criteria facet_criteria(30., 0.3, 0.03);
    Cell_criteria cell_criteria(3., 0.4);

    Criteria criteria(
        edge_criteria,
        facet_criteria,
        cell_criteria);

    C3t3 c3t3 =
        CGAL::make_mesh_3<C3t3>(
            domain,
            criteria,
            CGAL::parameters::no_exude().no_perturb());

    CGAL::remove_far_points_in_mesh_3(c3t3);

    assert(c3t3.number_of_facets_in_complex() > 0);
    assert(c3t3.number_of_edges_in_complex() > 0);

    Projector projector(domain);

    // Check one real Mesh_3 curve query. Again, only coarse geometric
    // consistency is checked.
    {
        const auto edge = *c3t3.edges_in_complex().begin();
        const Segment_3 segment =
            c3t3.triangulation().segment(edge);

        const std::array<Point_3, 2> points = {
            segment.source(),
            segment.target()
        };

        const auto tangent_space =
            projector.curve_edge_projection_line(
                std::make_pair(
                    c3t3.curve_index(edge),
                    edge),
                points);

        assert(
            tangent_space.projection_mode() !=
            Projection_mode::NONE);

        assert(finite_point(tangent_space.origin()));
        assert(tangent_space.vector().squared_length() > FT(0));

        assert(
            distance_to_cube_edges(tangent_space.origin()) <
            CURVE_TOLERANCE);
    }

    const auto status =
        CGAL::boundary_aware_mesh_smoothing(
            c3t3,
            projector);

    assert_all_finite(c3t3);
    assert(status.valid_mesh());
    assert(status.nb_invalid_elements == 0);

    assert(
        max_surface_distance(c3t3, domain) <
        SURFACE_TOLERANCE);

    assert(
        max_curve_distance_to_cube(c3t3) <
        CURVE_TOLERANCE);
}


// -----------------------------------------------------------------------------
// Signed_distance_function_projector
// -----------------------------------------------------------------------------

FT sphere_function(const Point_3& p)
{
    return
        p.x() * p.x() +
        p.y() * p.y() +
        p.z() * p.z() -
        FT(1);
}

std::pair<FT, Vector_3>
sphere_signed_distance(const Point_3& p)
{
    const FT squared_radius =
        p.x() * p.x() +
        p.y() * p.y() +
        p.z() * p.z();

    const FT radius = CGAL::sqrt(squared_radius);

    CGAL_precondition(radius != FT(0));

    return {
        radius - FT(1),
        Vector_3(
            p.x() / radius,
            p.y() / radius,
            p.z() / radius)
    };
}

template<typename C3t3>
double max_sphere_distance(const C3t3& c3t3)
{
    const auto& tr = c3t3.triangulation();

    double max_distance = 0.;

    for(const auto& facet : c3t3.facets_in_complex())
    {
        const Point_3 center =
            CGAL::centroid(tr.triangle(facet));

        const double distance = std::abs(
            CGAL::to_double(
                sphere_signed_distance(center).first));

        max_distance = (std::max)(
            max_distance,
            distance);
    }

    return max_distance;
}

void test_signed_distance_function_projector_iterations()
{
    std::size_t nb_evaluations = 0;

    auto sdf =
        [&nb_evaluations](const Point_3& p)
            -> std::pair<FT, Vector_3>
        {
            ++nb_evaluations;
            return sphere_signed_distance(p);
        };

    using Projector =
        CGAL::Mesh_smoothing_3::
            Signed_distance_function_projector<
                K,
                decltype(sdf)>;

    Projector projector(
        sdf,
        10,
        FT(1e-8));

    // Centroid = (2, 0, 0). For an exact SDF, the first iteration lands
    // directly on the sphere. The projector must still evaluate the field
    // again to confirm positional convergence.
    const std::vector<Point_3> points = {
        Point_3(2., -0.1, -0.1),
        Point_3(2.,  0.1, -0.1),
        Point_3(2.,  0.0,  0.2)
    };

    const auto tangent_space =
        projector.patch_face_projection_plane(0, points);

    assert(nb_evaluations >= 2);
    assert(finite_point(tangent_space.origin()));
    assert(tangent_space.vector().squared_length() > FT(0));

    // Deliberately loose: this is only a sanity check on the projection.
    assert(
        std::abs(
            CGAL::to_double(
                sphere_signed_distance(
                    tangent_space.origin()).first)) < 0.1);

    const std::array<Point_3, 2> curve_points = {
        Point_3(1., 0., 0.),
        Point_3(0., 1., 0.)
    };

    const auto curve =
        projector.curve_edge_projection_line(
            0,
            curve_points);

    assert(
        curve.projection_mode() ==
        Projection_mode::NONE);
}

void test_signed_distance_function_projector()
{
    using Domain = CGAL::Labeled_mesh_domain_3<K>;

    using Tr =
        typename CGAL::Mesh_triangulation_3<
            Domain,
            typename CGAL::Kernel_traits<Domain>::Kernel>::type;

    using C3t3 =
        CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

    using Criteria = CGAL::Mesh_criteria_3<Tr>;
    using Facet_criteria = typename Criteria::Facet_criteria;
    using Cell_criteria = typename Criteria::Cell_criteria;
    using Surface_patch_index =
        typename Domain::Surface_patch_index;

    namespace params = CGAL::parameters;

    Domain domain =
        Domain::create_implicit_mesh_domain(
            params::function = sphere_function,
            params::bounding_object =
                K::Sphere_3(CGAL::ORIGIN, FT(2)),
            params::p_rng = &CGAL::get_default_random(),
            params::relative_error_bound = 1e-3);

    Facet_criteria facet_criteria(
        0.,
        0.,
        0.3);

    Cell_criteria cell_criteria(
        0.,
        0.5);

    Criteria criteria(
        facet_criteria,
        cell_criteria);

    // Same robust initialization pattern as the Mesh_3 implicit-domain test.
    const std::vector<Point_3> initial_points = {
        Point_3( 1,  0,  0),
        Point_3( 0,  1,  0),
        Point_3( 0,  0,  1),
        Point_3(-1,  0,  0),
        Point_3( 0, -1,  0),
        Point_3( 0,  0, -1)
    };

    C3t3 c3t3;

    c3t3.insert_surface_points(
        initial_points.begin(),
        initial_points.end(),
        domain.index_from_surface_patch_index(
            Surface_patch_index(0, 1)));

    CGAL::refine_mesh_3(
        c3t3,
        domain,
        criteria,
        params::no_exude().no_perturb());

    CGAL::remove_far_points_in_mesh_3(c3t3);

    assert(c3t3.number_of_facets_in_complex() > 0);

    auto sdf =
        [](const Point_3& p)
            -> std::pair<FT, Vector_3>
        {
            return sphere_signed_distance(p);
        };

    using Projector =
        CGAL::Mesh_smoothing_3::
            Signed_distance_function_projector<
                K,
                decltype(sdf)>;

    Projector projector(
        sdf,
        10,
        FT(1e-8));

    const auto status =
        CGAL::boundary_aware_mesh_smoothing(
            c3t3,
            projector);

    assert_all_finite(c3t3);
    assert(status.valid_mesh());
    assert(status.nb_invalid_elements == 0);

    // Do not require exact projection: the target term is balanced against
    // tetrahedral quality during smoothing.
    assert(
        max_sphere_distance(c3t3) <
        IMPLICIT_TOLERANCE);
}

} // namespace

int main()
{
    test_polyhedral_mesh_domain_projector();
    test_polyhedral_mesh_domain_with_features_projector();

    test_signed_distance_function_projector_iterations();
    test_signed_distance_function_projector();

    return EXIT_SUCCESS;
}