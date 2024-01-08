// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef AOS_DEFS_H
#define AOS_DEFS_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Vector_3.h>

//#define USE_EPIC

#ifdef USE_EPIC
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
#else
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
#endif

using Geom_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
using Point = Geom_traits::Point_2;
using Curve = Geom_traits::Curve_2;
using Topol_traits = CGAL::Arr_spherical_topology_traits_2<Geom_traits>;
using Arrangement = CGAL::Arrangement_on_surface_2<Geom_traits, Topol_traits>;

// the following is from "arr_inexact_construction_segments.h":
using Segment = Geom_traits::X_monotone_curve_2;
using Vertex_handle = Arrangement::Vertex_handle;

// COUNTRIES AOS for grouping the faces by the country name
using Countries_dcel = CGAL::Arr_face_extended_dcel<Geom_traits, std::string>;
using Countries_topol_traits =
  CGAL::Arr_spherical_topology_traits_2<Geom_traits, Countries_dcel>;
using Countries_arr =
  CGAL::Arrangement_on_surface_2<Geom_traits, Countries_topol_traits>;


#endif
