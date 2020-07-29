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

#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

template <class K>
int test()
{
  // types
  typedef typename K::FT FT;
  typedef typename K::Line_3 Line;
  typedef typename K::Point_3 Point;
  typedef typename K::Segment_3 Segment;

  // load polyhedron
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  Polyhedron polyhedron;
  std::ifstream ifs("./data/tetrahedron.off");
  ifs >> polyhedron;

  // construct tree from facets
  typedef typename CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
  typedef typename CGAL::AABB_traits<K,Primitive> Traits;
  typedef typename CGAL::AABB_tree<Traits> Tree;
  typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;
  Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);

  // segment intersection query
  Point p((FT)-0.25,  (FT)0.251, (FT)0.255);
  Point q((FT) 0.25,  (FT)0.253, (FT)0.256);
  Segment pq(p,q);

  if(!tree.do_intersect(pq))
  {
    std::cerr << "no intersection found" << std::endl;
    return EXIT_FAILURE;
  }

  if(tree.number_of_intersected_primitives(pq) != 1)
  {
    std::cerr << "number of intersections different than one" << std::endl;
    return EXIT_FAILURE;
  }

  boost::optional<Object_and_primitive_id> any;
  any = tree.any_intersection(pq);
  if(!any)
  {
    std::cerr << "did not find any intersection" << std::endl;
    return EXIT_FAILURE;
  }
  Object_and_primitive_id op = *any;
  CGAL::Object object = op.first;
  Point point;
  if(CGAL::assign(point,object))
  {
    std::cout << "Intersection point: " << point << std::endl;
  }
  else
  {
    std::cerr << "intersection does not assign to a point" << std::endl;
    return EXIT_FAILURE;
  }

  // line intersection query
  Line line_pq(p,q);
  if(!tree.do_intersect(line_pq))
  {
    std::cerr << "no intersection found with line" << std::endl;
    return EXIT_FAILURE;
  }
  if(tree.number_of_intersected_primitives(line_pq) != 2)
  {
    std::cerr << "number of intersections different than two with line" << std::endl;
    return EXIT_FAILURE;
  }

  // closest point query
  Point r((FT)0.0, (FT)0.0, (FT)3.0);
  Point closest((FT)0.0, (FT)0.0, (FT)1.0);
  Point result = tree.closest_point(r);
  if(result != closest)
  {
    std::cerr << "wrong closest point" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int main()
{
  if(test<CGAL::Simple_cartesian<float> >() == EXIT_FAILURE)
    return EXIT_FAILURE;

  if(test<CGAL::Simple_cartesian<double> >() == EXIT_FAILURE)
    return EXIT_FAILURE;

  if(test<CGAL::Exact_predicates_inexact_constructions_kernel>() == EXIT_FAILURE)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}

/***EMACS SETTINGS***/
/* Local Variables: */
/* tab-width: 2     */
/* End:             */

