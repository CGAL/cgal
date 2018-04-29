// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Maxime Gimeno

// List of named parameters that we use in CGAL
CGAL_add_named_parameter(vertex_point_t, vertex_point, vertex_point_map)
CGAL_add_named_parameter(halfedge_index_t, halfedge_index, halfedge_index_map)
CGAL_add_named_parameter(edge_index_t, edge_index, edge_index_map)
CGAL_add_named_parameter(face_index_t, face_index, face_index_map)

CGAL_add_named_parameter(edge_is_constrained_t, edge_is_constrained, edge_is_constrained_map)
CGAL_add_named_parameter(first_index_t, first_index, first_index)
CGAL_add_named_parameter(number_of_iterations_t, number_of_iterations, number_of_iterations)
CGAL_add_named_parameter(verbosity_level_t, verbosity_level, verbosity_level)

CGAL_add_named_parameter(metis_options_t, METIS_options, METIS_options)
CGAL_add_named_parameter(vertex_partition_id_t, vertex_partition_id, vertex_partition_id_map)
CGAL_add_named_parameter(face_partition_id_t, face_partition_id, face_partition_id_map)

// List of named parameters that we use in the package 'Mesh_3'
CGAL_add_named_parameter(vertex_feature_degree_t, vertex_feature_degree, vertex_feature_degree_map)

// List of named parameters used in the package 'Polygon Mesh Processing'
CGAL_add_named_parameter(geom_traits_t, geom_traits, geom_traits)
CGAL_add_named_parameter(vertex_incident_patches_t, vertex_incident_patches, vertex_incident_patches_map)
CGAL_add_named_parameter(density_control_factor_t, density_control_factor, density_control_factor)
CGAL_add_named_parameter(use_delaunay_triangulation_t, use_delaunay_triangulation, use_delaunay_triangulation)
CGAL_add_named_parameter(fairing_continuity_t, fairing_continuity, fairing_continuity)
CGAL_add_named_parameter(sparse_linear_solver_t, sparse_linear_solver, sparse_linear_solver)
CGAL_add_named_parameter(number_of_relaxation_steps_t, number_of_relaxation_steps, number_of_relaxation_steps)
CGAL_add_named_parameter(protect_constraints_t, protect_constraints, protect_constraints)
CGAL_add_named_parameter(relax_constraints_t, relax_constraints, relax_constraints)
CGAL_add_named_parameter(vertex_is_constrained_t, vertex_is_constrained, vertex_is_constrained_map)
CGAL_add_named_parameter(face_patch_t, face_patch, face_patch_map)
CGAL_add_named_parameter(random_uniform_sampling_t, random_uniform_sampling, use_random_uniform_sampling)
CGAL_add_named_parameter(grid_sampling_t, grid_sampling, use_grid_sampling)
CGAL_add_named_parameter(monte_carlo_sampling_t, monte_carlo_sampling, use_monte_carlo_sampling)
CGAL_add_named_parameter(do_sample_edges_t, do_sample_edges, do_sample_edges)
CGAL_add_named_parameter(do_sample_vertices_t, do_sample_vertices, do_sample_vertices)
CGAL_add_named_parameter(do_sample_faces_t, do_sample_faces, do_sample_faces)
CGAL_add_named_parameter(number_of_points_on_faces_t, number_of_points_on_faces, number_of_points_on_faces)
CGAL_add_named_parameter(number_of_points_per_face_t, number_of_points_per_face, number_of_points_per_face)
CGAL_add_named_parameter(grid_spacing_t, grid_spacing, grid_spacing)
CGAL_add_named_parameter(number_of_points_per_edge_t, number_of_points_per_edge, number_of_points_per_edge)
CGAL_add_named_parameter(number_of_points_on_edges_t, number_of_points_on_edges, number_of_points_on_edges)
CGAL_add_named_parameter(nb_points_per_area_unit_t, nb_points_per_area_unit, number_of_points_per_area_unit)
CGAL_add_named_parameter(nb_points_per_distance_unit_t, nb_points_per_distance_unit, number_of_points_per_distance_unit)
CGAL_add_named_parameter(outward_orientation_t, outward_orientation, outward_orientation)
CGAL_add_named_parameter(overlap_test_t, overlap_test, do_overlap_test_of_bounded_sides)
CGAL_add_named_parameter(preserve_genus_t, preserve_genus, preserve_genus)

// List of named parameters that we use in the package 'Surface Mesh Simplification'
CGAL_add_named_parameter(get_cost_policy_t, get_cost_policy, get_cost)
CGAL_add_named_parameter(get_placement_policy_t, get_placement_policy, get_placement)

//to be documented
CGAL_add_named_parameter(face_normal_t, face_normal, face_normal_map)
CGAL_add_named_parameter(random_seed_t, random_seed, random_seed)
CGAL_add_named_parameter(do_project_t, do_project, do_project)

//internal
CGAL_add_named_parameter(weight_calculator_t, weight_calculator, weight_calculator)

// List of named parameters used in the Point Set Processing package
CGAL_add_named_parameter(point_t, point_map, point_map)
CGAL_add_named_parameter(query_point_t, query_point_map, query_point_map)
CGAL_add_named_parameter(normal_t, normal_map, normal_map)
CGAL_add_named_parameter(diagonalize_traits_t, diagonalize_traits, diagonalize_traits)
CGAL_add_named_parameter(svd_traits_t, svd_traits, svd_traits)
CGAL_add_named_parameter(sharpness_angle_t, sharpness_angle, sharpness_angle)
CGAL_add_named_parameter(edge_sensitivity_t, edge_sensitivity, edge_sensitivity)
CGAL_add_named_parameter(neighbor_radius_t, neighbor_radius, neighbor_radius)
CGAL_add_named_parameter(number_of_output_points_t, number_of_output_points, number_of_output_points)
CGAL_add_named_parameter(size_t, size, size)
CGAL_add_named_parameter(maximum_variation_t, maximum_variation, maximum_variation)
CGAL_add_named_parameter(degree_fitting_t, degree_fitting, degree_fitting)
CGAL_add_named_parameter(degree_monge_t, degree_monge, degree_monge)
CGAL_add_named_parameter(threshold_percent_t, threshold_percent, threshold_percent)
CGAL_add_named_parameter(threshold_distance_t, threshold_distance, threshold_distance)
CGAL_add_named_parameter(attraction_factor_t, attraction_factor, attraction_factor)
CGAL_add_named_parameter(plane_t, plane_map, plane_map)
CGAL_add_named_parameter(plane_index_t, plane_index_map, plane_index_map)
CGAL_add_named_parameter(select_percentage_t, select_percentage, select_percentage)
CGAL_add_named_parameter(require_uniform_sampling_t, require_uniform_sampling, require_uniform_sampling)
