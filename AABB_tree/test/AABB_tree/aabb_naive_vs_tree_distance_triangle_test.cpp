// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <fstream>
#include <iostream>

#include <CGAL/Timer.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>

#include "AABB_test_util.h"


template<class K, class Tree, class Polyhedron, Primitive_type Type>
void test_impl(Tree& tree, Polyhedron& p, const double duration)
{
  tree.accelerate_distance_queries(p.points_begin(),p.points_end());

  typedef Tree_vs_naive<Tree, Polyhedron, K, Type> Tester;
  Tester tester(tree, p);

  tester.test_all_distance_methods(duration);
}

int main(void)
{
  std::cout << "AABB naive vs tree distance (triangle primitive) tests" << std::endl;

  const double duration = 0.1;
  test_kernels<TRIANGLE>("data/cube.off",duration);
  test_kernels<TRIANGLE>("data/coverrear.off",duration);
  test_kernels<TRIANGLE>("data/finger.off",duration);
  test_kernels<TRIANGLE>(CGAL::data_file_path("meshes/pinion_small.off"),duration);

  return EXIT_SUCCESS;
}
