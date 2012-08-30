// Copyright (c) 2011 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://scm.gforge.inria.fr/svn/cgal/branches/features/Mesh_3-experimental-GF/Mesh_3/include/CGAL/Mesh_3/config.h $
// $Id: config.h 70893 2012-07-31 10:43:20Z jtournoi $
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_CONFIG_H
#define CGAL_MESH_3_CONFIG_H

#include <CGAL/Mesh_3/config.h>

#include <CGAL/config.h>


// Use optimisations of Mesh_3
#  define CGAL_COMPACT_MESH_VERTEX_CELL
#  define CGAL_FREEZE_VERTICES
#  define CGAL_IMPROVE_FREEZE
#  define CGAL_REGULAR_TRIANGULATION_3_USE_CIRCUMCENTER_CACHE
#  define CGAL_INTRUSIVE_LIST
#  define CGAL_CONSTRUCT_INTRUSIVE_LIST_RANGE_CONSTRUCTOR
#  define CGAL_MESH_3_NEW_GET_FACETS
#  define CGAL_MESH_3_GET_FACETS_USING_INTRUSIVE_LIST
#  define CGAL_MESH_3_SIZING_FIELD_INEXACT_LOCATE
#  define FORCE_STRUCTURAL_FILTERING

//experimental
#  define CGAL_FASTER_BUILD_QUEUE

// to disable new features :
//#  define CGAL_TDS_3_NOT_JUST_INCIDENT_CELLS

//should not be used
//#define CGAL_MESH_3_OLD_MINIMUM_DIHEDRAL_ANGLE

//experimental
//#define CGAL_MESH_3_PROTECTION_NON_LINEAR
//#define CGAL_POLYHEDRAL_MESH_DOMAIN_USE_GRID
//#define CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS


// CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS
// implies CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
#ifdef CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS
#  ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
#    define CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3 1
#  endif
#endif // CGAL_MESH_3_NEW_ROBUST_INTERSECTION_TRAITS

// CGAL_MESH_3_VERBOSE implies CGAL_MESH_3_OPTIMIZER_VERBOSE
#ifdef CGAL_MESH_3_VERBOSE
#  ifndef CGAL_MESH_3_OPTIMIZER_VERBOSE
#    define CGAL_MESH_3_OPTIMIZER_VERBOSE 1
#  endif
#endif

#endif // CGAL_MESH_3_CONFIG_H
