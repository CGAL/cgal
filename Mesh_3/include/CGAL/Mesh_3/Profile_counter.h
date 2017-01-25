// Copyright (c) 2005,2006,2008  INRIA Sophia-Antipolis (France).
// Copyright (c) 2011 GeometryFactory Sarl (France)
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
// $URL: https://scm.gforge.inria.fr/svn/cgal/branches/features/Mesh_3-experimental-GF/Mesh_3/include/CGAL/Mesh_3/Profile_counter.h $
// $Id: Profile_counter.h 66846 2011-12-16 13:48:50Z lrineau $
//
//
// Author(s)     : Sylvain Pion, Laurent Rineau

#ifndef CGAL_MESH_3_PROFILE_COUNTER_H
#define CGAL_MESH_3_PROFILE_COUNTER_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/Profile_counter.h>

#ifdef CGAL_MESH_3_PROFILE
#  define CGAL_MESH_3_PROFILER(Y) \
          { static CGAL::Profile_counter tmp(Y); ++tmp; }
#  define CGAL_MESH_3_HISTOGRAM_PROFILER(Y, Z) \
          { static CGAL::Profile_histogram_counter tmp(Y); tmp(Z); }
#  define CGAL_MESH_3_BRANCH_PROFILER(Y, NAME) \
          static CGAL::Profile_branch_counter NAME(Y); ++NAME;
#  define CGAL_MESH_3_BRANCH_PROFILER_BRANCH(NAME) \
          NAME.increment_branch();
#  define CGAL_MESH_3_BRANCH_PROFILER_3(Y, NAME) \
          static CGAL::Profile_branch_counter_3 NAME(Y); ++NAME;
#  define CGAL_MESH_3_BRANCH_PROFILER_BRANCH_1(NAME) \
          NAME.increment_branch_1();
#  define CGAL_MESH_3_BRANCH_PROFILER_BRANCH_2(NAME) \
          NAME.increment_branch_2();
#else
#  define CGAL_MESH_3_PROFILER(Y)
#  define CGAL_MESH_3_HISTOGRAM_PROFILER(Y, Z)
#  define CGAL_MESH_3_BRANCH_PROFILER(Y, NAME)
#  define CGAL_MESH_3_BRANCH_PROFILER_BRANCH(NAME)
#  define CGAL_MESH_3_BRANCH_PROFILER_3(Y, NAME) 
#  define CGAL_MESH_3_BRANCH_PROFILER_BRANCH_1(NAME)
#  define CGAL_MESH_3_BRANCH_PROFILER_BRANCH_2(NAME) 
#endif

#endif // CGAL_MESH_3_PROFILE_COUNTER_H
