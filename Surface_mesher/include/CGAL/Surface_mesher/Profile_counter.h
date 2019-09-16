// Copyright (c) 2008  GeometryFactory (France)
// All rights reserved.
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
// Author(s)     : Laurent Rineau

// This file is an adaptation of <CGAL/Profile_timer.h>, so that the
// macros are prefixed with CGAL_SURFACE_MESHER_ instead of CGAL_.

#ifndef CGAL_SURFACE_MESHER_PROFILE_COUNTER_H
#define CGAL_SURFACE_MESHER_PROFILE_COUNTER_H

#include <CGAL/Profile_counter.h>

#ifdef CGAL_SURFACE_MESHER_PROFILE
#  define CGAL_SURFACE_MESHER_PROFILER(Y) \
          { static CGAL::Profile_counter tmp(Y); ++tmp; }
#  define CGAL_SURFACE_MESHER_HISTOGRAM_PROFILER(Y, Z) \
          { static CGAL::Profile_histogram_counter tmp(Y); tmp(Z); }
#  define CGAL_SURFACE_MESHER_BRANCH_PROFILER(Y, NAME) \
          static CGAL::Profile_branch_counter NAME(Y); ++NAME;
#  define CGAL_SURFACE_MESHER_BRANCH_PROFILER_BRANCH(NAME) \
          NAME.increment_branch();
#  define CGAL_SURFACE_MESHER_BRANCH_PROFILER_3(Y, NAME) \
          static CGAL::Profile_branch_counter_3 NAME(Y); ++NAME;
#  define CGAL_SURFACE_MESHER_BRANCH_PROFILER_BRANCH_1(NAME) \
          NAME.increment_branch_1();
#  define CGAL_SURFACE_MESHER_BRANCH_PROFILER_BRANCH_2(NAME) \
          NAME.increment_branch_2();
#else
#  define CGAL_SURFACE_MESHER_PROFILER(Y)
#  define CGAL_SURFACE_MESHER_HISTOGRAM_PROFILER(Y, Z)
#  define CGAL_SURFACE_MESHER_BRANCH_PROFILER(Y, NAME)
#  define CGAL_SURFACE_MESHER_BRANCH_PROFILER_BRANCH(NAME)
#  define CGAL_SURFACE_MESHER_BRANCH_PROFILER_3(Y, NAME) 
#  define CGAL_SURFACE_MESHER_BRANCH_PROFILER_BRANCH_1(NAME)
#  define CGAL_SURFACE_MESHER_BRANCH_PROFILER_BRANCH_2(NAME) 
#endif

#endif // CGAL_SURFACE_MESHER_PROFILE_COUNTER_H
