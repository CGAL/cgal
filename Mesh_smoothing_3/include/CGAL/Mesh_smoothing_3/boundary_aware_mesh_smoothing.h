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
#include <CGAL/Mesh_smoothing_3/Smoothing_status.h>
#include <CGAL/Mesh_smoothing_3/projectors.h>

#include <CGAL/Named_function_parameters.h>

#include <array>
#include <vector>
#include <tuple>

namespace CGAL {

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
*                           Algorithm will stop before if it reaches convergence.
*                           Untangling usually requires more iterations (up to thousands) for hard cases. }
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`100`}
*   \cgalParamNEnd
*   \cgalParamNBegin{verbose}
*     \cgalParamDescription{If `true`, the function prints information about the smoothing
*                           process to the standard output.}
*     \cgalParamType{`bool`}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
*   \cgalParamNBegin{maximum_running_time}
*     \cgalParamDescription{Maximum allowed time for the smoothing process in seconds.}
*     \cgalParamType{`double`}
*     \cgalParamDefault{`-1`}
*     \cgalParamExtra{Pre-processing will not be stopped and negative times will be ignored.}
*   \cgalParamNEnd
*   \cgalParamNBegin{max_number_of_evaluations}
*     \cgalParamDescription{Maximum number of quality metric evaluations for smoothing.}
*     \cgalParamType{`int`}
*     \cgalParamDefault{`-1`}
*     \cgalParamExtra{Strongly correlated to running time but will scale linearly with mesh size.}
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
*
*
* \return a `Mesh_smoothing_3::Smoothing_status` object containing information about the smoothing process and convergence.
*
*/
template<typename C3t3,
         typename CTS,
         typename NamedParameters = parameters::Default_named_parameters>
Mesh_smoothing_3::Smoothing_status boundary_aware_mesh_smoothing  (
    C3t3& c3t3,
    CTS const & cts,
    NamedParameters const & np = parameters::default_values()
)
{
    Mesh_smoothing_3::Smoothing_status return_status;
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


    int max_nb_metric_evaluations = choose_parameter(get_parameter(np, internal_np::max_number_of_evaluations), -1);

    double time_limit = choose_parameter(get_parameter(np, internal_np::maximum_running_time), -1.);

    // C3t3_smoother will automatically mark facets and features
    CGAL::Mesh_smoothing_3::C3t3_smoother smoother(c3t3);

    // locks through property maps
    for (auto v : c3t3.triangulation().finite_vertex_handles()) {
        if (!get(vcmap, v)) continue;
        smoother.set_vertex_lock(v, true);
    }

    for (const auto& e : c3t3.triangulation().finite_edges()) {
        const auto vertices = c3t3.triangulation().vertices(e);
        const std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle> edge_01{vertices[0], vertices[1]};
        const std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle> edge_10{vertices[1], vertices[0]};
        if (get(ecmap, edge_01) || get(ecmap, edge_10)) {
            smoother.set_vertex_lock(vertices[0], true);
            smoother.set_vertex_lock(vertices[1], true);
        }
    }

    for (const auto& f : c3t3.triangulation().finite_facets()) {
       const auto mirror_f = c3t3.triangulation().mirror_facet(f);
       if (get(fcmap, f) || get(fcmap, mirror_f)) {
           for (const auto& v : c3t3.triangulation().vertices(f))
               smoother.set_vertex_lock(v, true);
       }
    }

    // locks corners
    for (auto c : c3t3.vertices_in_complex()) {
        smoother.set_vertex_lock(c, true);
    }

    using Point_3 = typename  C3t3::Triangulation::Geom_traits::Point_3;
    using Tangent_space = typename CTS::Tangent_space;

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

    smoother.set_predicates_mode(Mesh_smoothing_3::Parameters::STRONG_ENFORCEMENT);
    smoother.set_verbose(verbose);
    smoother.set_max_number_of_iteration(max_iterations);
    smoother.set_maximum_running_time(time_limit);
    smoother.set_maximum_number_of_metric_evaluations(max_nb_metric_evaluations);

    return_status.add_time(true);

    smoother.set_input_smoothing_status(return_status);

    return_status = smoother.run();

    return return_status;
}

}

#endif // CGAL_MESH_SMOOTHING_3_BOUNDARY_AWARE_MESH_SMOOTHING_H
