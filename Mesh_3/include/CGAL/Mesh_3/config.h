// Copyright (c) 2011 GeometryFactory Sarl (France).
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
// $URL: https://scm.gforge.inria.fr/svn/cgal/branches/features/Mesh_3-experimental-GF/Mesh_3/include/CGAL/Mesh_3/config.h $
// $Id: config.h 70893 2012-07-31 10:43:20Z jtournoi $
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_CONFIG_H
#define CGAL_MESH_3_CONFIG_H 1

#include <CGAL/Mesh_3/config.h>

#include <CGAL/config.h>

//#define CGAL_MESH_3_VERBOSE 1

// Use optimisations of Mesh_3
#  define CGAL_REGULAR_TRIANGULATION_3_USE_CIRCUMCENTER_CACHE 1
#  define CGAL_INTRUSIVE_LIST 1
#  define CGAL_CONSTRUCT_INTRUSIVE_LIST_RANGE_CONSTRUCTOR 1
#  define CGAL_MESH_3_NEW_GET_FACETS 1
#  define CGAL_MESH_3_GET_FACETS_USING_INTRUSIVE_LIST 1
#  define CGAL_MESH_3_SIZING_FIELD_INEXACT_LOCATE 1
#  define FORCE_STRUCTURAL_FILTERING 1
#  define CGAL_NEW_INCIDENT_SLIVERS 1

//experimental
#  define CGAL_FASTER_BUILD_QUEUE 1

//should not be used
//#define CGAL_MESH_3_OLD_MINIMUM_DIHEDRAL_ANGLE 1

//experimental
//#define CGAL_MESH_3_PROTECTION_NON_LINEAR 1
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
