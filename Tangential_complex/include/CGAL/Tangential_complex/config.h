// Copyright (c) 2014  INRIA Sophia-Antipolis (France)
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Clement Jamin

#ifndef CGAL_TC_CONFIG_H
#define CGAL_TC_CONFIG_H

#include <CGAL/config.h>

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false
// positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

//=========================== Alpha-TC or not? ================================

//#define CGAL_ALPHA_TC
  //#define CGAL_USE_A_FIXED_ALPHA
  const double CGAL_TC_ALPHA_VALUE = 0.3;

//========================= Debugging & profiling =============================
#define CGAL_TC_PROFILING
#define CGAL_TC_VERBOSE
//#define CGAL_TC_VERY_VERBOSE
//#define CGAL_TC_PERFORM_EXTRA_CHECKS
#define CGAL_TC_SHOW_DETAILED_STATS_FOR_INCONSISTENCIES

// Solving inconsistencies: only perturb the vertex, the simplex or more?
//#define CGAL_TC_PERTURB_THE_CENTER_VERTEX_ONLY
//#define CGAL_TC_PERTURB_THE_SIMPLEX_ONLY
//#define CGAL_TC_PERTURB_THE_1_STAR
#define CGAL_TC_PERTURB_N_CLOSEST_POINTS // perturb the CGAL_TC_NUMBER_OF_PERTURBED_POINTS closest points
// Otherwise, perturb one random point of the simplex

// Only used if CGAL_TC_PERTURB_N_CLOSEST_POINTS is defined
#define CGAL_TC_NUMBER_OF_PERTURBED_POINTS(intr_dim) (1)
//#define CGAL_TC_NUMBER_OF_PERTURBED_POINTS(intr_dim) (intr_dim + 2)

#define CGAL_MESH_D_PROFILING
#define CGAL_MESH_D_VERBOSE

//========================= Strategy ==========================================
#define CGAL_TC_NANOFLANN_IS_AVAILABLE
//#define CGAL_TC_ANN_IS_AVAILABLE
//#define CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER
#define CGAL_TC_GLOBAL_REFRESH
//#define CGAL_TC_ON_DEMAND_REFRESH // CJTODO: not implemented yet
    // The idea is to perform a global refresh + some local refreshes, just
    // for local tri where there are some inconsistencies
    // But be careful: refreshing the TC may invalidate cells, so the
    // incident cells have to be recomputed again
#define CGAL_TC_PERTURB_POSITION
# define CGAL_TC_PERTURB_POSITION_TANGENTIAL // default
//# define CGAL_TC_PERTURB_POSITION_GLOBAL
//#define CGAL_TC_PERTURB_WEIGHT
//#define CGAL_TC_PERTURB_TANGENT_SPACE

//========================= Parameters ========================================

// PCA will use BASE_VALUE_FOR_PCA^intrinsic_dim points
const std::size_t BASE_VALUE_FOR_PCA = 5;
const unsigned int BASE_VALUE_FOR_ALPHA_TC_NEIGHBORHOOD = 5;

#endif // CGAL_TC_CONFIG_H
