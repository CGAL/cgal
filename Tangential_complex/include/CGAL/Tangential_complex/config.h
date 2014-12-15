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

//========================= Debugging & profiling =============================
#define CGAL_TC_PROFILING
#define CGAL_TC_VERBOSE

// Solving inconsistencies: only change the weights of the inconsistent simplex
// or more?
//#define CGAL_TC_ONLY_CHANGE_SIMPLEX_WEIGHTS

// Otherwise...
const int CGAL_TC_NUMBER_OF_ADDITIONNAL_PERTURBED_POINTS = 1;

//========================= Strategy ==========================================
//#define CGAL_TC_USE_NANOFLANN
//#define CGAL_TC_USE_SLOW_BUT_ACCURATE_SPARSIFIER

//========================= Parameters ========================================
const std::size_t NUM_POINTS_FOR_PCA = 50;
const double INPUT_SPARSITY = 0.1;

#endif // CGAL_TC_CONFIG_H
