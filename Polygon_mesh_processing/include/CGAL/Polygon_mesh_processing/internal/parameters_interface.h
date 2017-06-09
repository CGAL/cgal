// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Maxime Gimeno

// List of named parameters used in the Polygon Mesh Processing package

CGAL_add_named_parameter(geom_traits_t, geom_traits, geom_traits)
CGAL_add_named_parameter(density_control_factor_t, density_control_factor, density_control_factor)
CGAL_add_named_parameter(use_delaunay_triangulation_t, use_delaunay_triangulation, use_delaunay_triangulation)
CGAL_add_named_parameter(fairing_continuity_t, fairing_continuity, fairing_continuity)
CGAL_add_named_parameter(sparse_linear_solver_t, sparse_linear_solver, sparse_linear_solver)
CGAL_add_named_parameter(number_of_iterations_t, number_of_iterations, number_of_iterations)
CGAL_add_named_parameter(number_of_relaxation_steps_t, number_of_relaxation_steps, number_of_relaxation_steps)
CGAL_add_named_parameter(protect_constraints_t, protect_constraints, protect_constraints)
CGAL_add_named_parameter(relax_constraints_t, relax_constraints, relax_constraints)
CGAL_add_named_parameter(vertex_is_constrained_t, vertex_is_constrained, vertex_is_constrained_map)
CGAL_add_named_parameter(face_patch_t, face_patch, face_patch_map)
CGAL_add_named_parameter(random_uniform_sampling_t, random_uniform_sampling, random_uniform_sampling)
CGAL_add_named_parameter(grid_sampling_t, grid_sampling, grid_sampling)
CGAL_add_named_parameter(monte_carlo_sampling_t, monte_carlo_sampling, monte_carlo_sampling)
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

//to be documented
CGAL_add_named_parameter(face_normal_t, face_normal, face_normal_map)
CGAL_add_named_parameter(random_seed_t, random_seed, random_seed)
CGAL_add_named_parameter(do_project_t, do_project, do_project)

//internal
CGAL_add_named_parameter(weight_calculator_t, weight_calculator, weight_calculator)
