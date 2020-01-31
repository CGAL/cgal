// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>


#include "AABB_test_util.h"


template<class K, class Tree, class Polyhedron, Primitive_type Type>
void test_impl(Tree& tree, Polyhedron& p, const double duration)
{
  tree.accelerate_distance_queries(p.points_begin(),p.points_end());

  test_distance_speed<Tree,K>(tree,duration);
  test_all_distance_query_types<Tree,K>(tree);
}

int main(void)
{
  std::cout << "AABB distance to edge tests" << std::endl;
  const double duration = 0.1;
  test_kernels<SEGMENT>("./data/cube.off",duration);
  test_kernels<SEGMENT>("./data/finger.off",duration);
  test_kernels<SEGMENT>("./data/pinion.off",duration);
  test_kernels<SEGMENT>("./data/coverrear.off",duration);
  return 0;
}
