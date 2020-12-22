// Copyright (c) 2011 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_CONFIG_H
#define CGAL_MESH_3_CONFIG_H 1

#include <CGAL/license/Triangulation_3.h>


#include <CGAL/config.h>

//#define CGAL_MESH_3_VERBOSE 1

// Use optimisations of Mesh_3
#  define CGAL_INTRUSIVE_LIST 1
#  define CGAL_CONSTRUCT_INTRUSIVE_LIST_RANGE_CONSTRUCTOR 1
#  define CGAL_MESH_3_NEW_GET_FACETS 1
#  define CGAL_MESH_3_GET_FACETS_USING_INTRUSIVE_LIST 1
#  define CGAL_MESH_3_SIZING_FIELD_INEXACT_LOCATE 1
#  define FORCE_STRUCTURAL_FILTERING 1
#  define CGAL_NEW_INCIDENT_SLIVERS 1

//experimental
#  define CGAL_FASTER_BUILD_QUEUE 1
//#  define CGAL_SEQUENTIAL_MESH_3_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE
//#  define CGAL_PARALLEL_MESH_3_DO_NOT_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE // slower / not recommended

//should not be used
//#define CGAL_MESH_3_OLD_MINIMUM_DIHEDRAL_ANGLE 1

//experimental
#define CGAL_MESH_3_NO_PROTECTION_NON_LINEAR 1
#define CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS 1


// CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS
// implies CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
#ifdef CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS
#  ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
#    ifndef CGAL_MESH_3_DEACTIVATE_NO_LONGER_CALLS_DO_INTERSECT_3
#      define CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3 1
#    endif
#  endif
#endif // CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS

// CGAL_MESH_3_VERBOSE implies CGAL_MESH_3_OPTIMIZER_VERBOSE
#ifdef CGAL_MESH_3_VERBOSE
#  ifndef CGAL_MESH_3_OPTIMIZER_VERBOSE
#    define CGAL_MESH_3_OPTIMIZER_VERBOSE 1
#  endif
#endif

#endif // CGAL_MESH_3_CONFIG_H
