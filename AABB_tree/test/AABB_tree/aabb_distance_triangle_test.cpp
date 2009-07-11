// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <iostream>
#include <fstream>
#include <CGAL/Timer.h>

#include <CGAL/AABB_intersections.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

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
  std::cout << "AABB distance to triangle tests" << std::endl;
  const double duration = 0.1;
  test_kernels<TRIANGLE>("./data/cube.off",duration);
  test_kernels<TRIANGLE>("./data/pinion.off",duration);
  test_kernels<TRIANGLE>("./data/finger.off",duration);
  test_kernels<TRIANGLE>("./data/coverrear.off",duration);
  return 0;
}
