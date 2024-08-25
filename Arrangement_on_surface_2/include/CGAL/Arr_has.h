// Copyright (c) 2005,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel      <efif@post.tau.ac.il>
//            Eric Berberich <ericb@post.tau.ac.il>

#ifndef CGAL_ARR_HAS_H
#define CGAL_ARR_HAS_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <type_traits>

namespace CGAL {

// Compare_x_2
// Helper trait to check for the presence of nested Compare_x_2
template <typename, typename = std::void_t<>>
struct has_compare_x_2 : std::false_type {};

// Specialization if the nested type Compare_x_2 exists
template <typename T>
struct has_compare_x_2<T, std::void_t<typename T::Compare_x_2>> : std::true_type {};

// Compare_xy_2
// Helper trait to check for the presence of nested Compare_xy_2
template <typename, typename = std::void_t<>>
struct has_compare_xy_2 : std::false_type {};

// Specialization if the nested type Compare_xy_2 exists
template <typename T>
struct has_compare_xy_2<T, std::void_t<typename T::Compare_xy_2>> : std::true_type {};

// Construct_min_vertex_2
// Helper trait to check for the presence of nested Construct_min_vertex_2
template <typename, typename = std::void_t<>>
struct has_construct_min_vertex_2 : std::false_type {};

// Specialization if the nested type Construct_min_vertex_2 exists
template <typename T>
struct has_construct_min_vertex_2<T, std::void_t<typename T::Construct_min_vertex_2>> : std::true_type {};

// Construct_max_vertex_2
// Helper trait to check for the presence of nested Construct_max_vertex_2
template <typename, typename = std::void_t<>>
struct has_construct_max_vertex_2 : std::false_type {};

// Specialization if the nested type Construct_max_vertex_2 exists
template <typename T>
struct has_construct_max_vertex_2<T, std::void_t<typename T::Construct_max_vertex_2>> : std::true_type {};

// Is_vertical_2
// Helper trait to check for the presence of nested Is_vertical_2
template <typename, typename = std::void_t<>>
struct has_is_vertical_2 : std::false_type {};

// Specialization if the nested type Is_vertical_2 exists
template <typename T>
struct has_is_vertical_2<T, std::void_t<typename T::Is_vertical_2>> : std::true_type {};

// Compare_y_at_x_2
// Helper trait to check for the presence of nested Compare_y_at_x_2
template <typename, typename = std::void_t<>>
struct has_compare_y_at_x_2 : std::false_type {};

// Specialization if the nested type Compare_y_at_x_2 exists
template <typename T>
struct has_compare_y_at_x_2<T, std::void_t<typename T::Compare_y_at_x_2>> : std::true_type {};

// Equal
// Helper trait to check for the presence of nested Equal
template <typename, typename = std::void_t<>>
struct has_equal_2 : std::false_type {};

// Specialization if the nested type Equal exists
template <typename T>
struct has_equal_2<T, std::void_t<typename T::Equal>> : std::true_type {};

// Compare_y_at_x_left_2
// Helper trait to check for the presence of nested Compare_y_at_x_left_2
template <typename, typename = std::void_t<>>
struct has_compare_y_at_x_left_2 : std::false_type {};

// Specialization if the nested type Compare_y_at_x_left_2 exists
template <typename T>
struct has_compare_y_at_x_left_2<T, std::void_t<typename T::Compare_y_at_x_left_2>> : std::true_type {};

// Compare_y_at_x_right_2
// Helper trait to check for the presence of nested Compare_y_at_x_right_2
template <typename, typename = std::void_t<>>
struct has_compare_y_at_x_right_2 : std::false_type {};

// Specialization if the nested type Compare_y_at_x_right_2 exists
template <typename T>
struct has_compare_y_at_x_right_2<T, std::void_t<typename T::Compare_y_at_x_right_2>> : std::true_type {};

// Make_x_monotone_2
// Helper trait to check for the presence of nested Make_x_monotone_2
template <typename, typename = std::void_t<>>
struct has_make_x_monotone_2 : std::false_type {};

// Specialization if the nested type Make_x_monotone_2 exists
template <typename T>
struct has_make_x_monotone_2<T, std::void_t<typename T::Make_x_monotone_2>> : std::true_type {};

// Split_2
// Helper trait to check for the presence of nested Split_2
template <typename, typename = std::void_t<>>
struct has_split_2 : std::false_type {};

// Specialization if the nested type Split_2 exists
template <typename T>
struct has_split_2<T, std::void_t<typename T::Split_2>> : std::true_type {};

// Intersect_2
// Helper trait to check for the presence of nested Intersect_2
template <typename, typename = std::void_t<>>
struct has_intersect_2 : std::false_type {};

// Specialization if the nested type Intersect_2 exists
template <typename T>
struct has_intersect_2<T, std::void_t<typename T::Intersect_2>> : std::true_type {};

// Are_mergeable_2
// Helper trait to check for the presence of nested Are_mergeable_2
template <typename, typename = std::void_t<>>
struct has_are_mergeable_2 : std::false_type {};

// Specialization if the nested type Are_mergeable_2 exists
template <typename T>
struct has_are_mergeable_2<T, std::void_t<typename T::Are_mergeable_2>> : std::true_type {};

// Merge_2
// Helper trait to check for the presence of nested Merge_2
template <typename, typename = std::void_t<>>
struct has_merge_2 : std::false_type {};

// Specialization if the nested type Merge_2 exists
template <typename T>
struct has_merge_2<T, std::void_t<typename T::Merge_2>> : std::true_type {};

// Construct_opposite_2
// Helper trait to check for the presence of nested Construct_opposite_2
template <typename, typename = std::void_t<>>
struct has_construct_opposite_2 : std::false_type {};

// Specialization if the nested type Construct_opposite_2 exists
template <typename T>
struct has_construct_opposite_2<T, std::void_t<typename T::Construct_opposite_2>> : std::true_type {};

// Compare_endpoints_xy_2
// Helper trait to check for the presence of nested Compare_endpoints_xy_2
template <typename, typename = std::void_t<>>
struct has_compare_endpoints_xy_2 : std::false_type {};

// Specialization if the nested type Compare_endpoints_xy_2 exists
template <typename T>
struct has_compare_endpoints_xy_2<T, std::void_t<typename T::Compare_endpoints_xy_2>> : std::true_type {};



// Approximate_2
// Helper trait to check for the presence of nested Approximate_2
template <typename, typename = std::void_t<>>
struct has_approximate_2 : std::false_type {};

// Specialization if the nested type Approximate_2 exists
template <typename T>
struct has_approximate_2<T, std::void_t<typename T::Approximate_2>> : std::true_type {};


// Parameter_space_in_x_2
// Helper trait to check for the presence of nested Parameter_space_in_x_2
template <typename, typename = std::void_t<>>
struct has_parameter_space_in_x_2 : std::false_type {};

// Specialization if the nested type Parameter_space_in_x_2 exists
template <typename T>
struct has_parameter_space_in_x_2<T, std::void_t<typename T::Parameter_space_in_x_2>> : std::true_type {};

// Is_on_x_identification_2
// Helper trait to check for the presence of nested Is_on_x_identification_2
template <typename, typename = std::void_t<>>
struct has_is_on_x_identification_2 : std::false_type {};

// Specialization if the nested type Is_on_x_identification_2 exists
template <typename T>
struct has_is_on_x_identification_2<T, std::void_t<typename T::Is_on_x_identification_2>> : std::true_type {};

// Compare_y_on_boundary_2
// Helper trait to check for the presence of nested Compare_y_on_boundary_2
template <typename, typename = std::void_t<>>
struct has_compare_y_on_boundary_2 : std::false_type {};

// Specialization if the nested type Compare_y_on_boundary_2 exists
template <typename T>
struct has_compare_y_on_boundary_2<T, std::void_t<typename T::Compare_y_on_boundary_2>> : std::true_type {};

// Compare_y_near_boundary_2
// Helper trait to check for the presence of nested Compare_y_near_boundary_2
template <typename, typename = std::void_t<>>
struct has_compare_y_near_boundary_2 : std::false_type {};

// Specialization if the nested type Compare_y_near_boundary_2 exists
template <typename T>
struct has_compare_y_near_boundary_2<T, std::void_t<typename T::Compare_y_near_boundary_2>> : std::true_type {};

// Parameter_space_in_y_2
// Helper trait to check for the presence of nested Parameter_space_in_y_2
template <typename, typename = std::void_t<>>
struct has_parameter_space_in_y_2 : std::false_type {};

// Specialization if the nested type Parameter_space_in_y_2 exists
template <typename T>
struct has_parameter_space_in_y_2<T, std::void_t<typename T::Parameter_space_in_y_2>> : std::true_type {};

// Is_on_y_identification_2
// Helper trait to check for the presence of nested Is_on_y_identification_2
template <typename, typename = std::void_t<>>
struct has_is_on_y_identification_2 : std::false_type {};

// Specialization if the nested type Is_on_y_identification_2 exists
template <typename T>
struct has_is_on_y_identification_2<T, std::void_t<typename T::Is_on_y_identification_2>> : std::true_type {};

// Compare_x_on_boundary_2
// Helper trait to check for the presence of nested Compare_x_on_boundary_2
template <typename, typename = std::void_t<>>
struct has_compare_x_on_boundary_2 : std::false_type {};

// Specialization if the nested type Compare_x_on_boundary_2 exists
template <typename T>
struct has_compare_x_on_boundary_2<T, std::void_t<typename T::Compare_x_on_boundary_2>> : std::true_type {};

// Compare_x_near_boundary_2
// Helper trait to check for the presence of nested Compare_x_near_boundary_2
template <typename, typename = std::void_t<>>
struct has_compare_x_near_boundary_2 : std::false_type {};

// Specialization if the nested type Compare_x_near_boundary_2 exists
template <typename T>
struct has_compare_x_near_boundary_2<T, std::void_t<typename T::Compare_x_near_boundary_2>> : std::true_type {};

}

#endif
