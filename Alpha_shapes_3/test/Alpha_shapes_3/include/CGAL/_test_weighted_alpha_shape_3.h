// Copyright (c) 1998-2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// $Date$
//
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_TEST_WEIGHTED_ALPHA_SHAPE_3_H
#define CGAL_TEST_WEIGHTED_ALPHA_SHAPE_3_H

#include <list>
#include <fstream>
#include "_count_alpha.h"

template <class Weighted_point>
bool
file_input(std::ifstream& is, std::list<Weighted_point>& L)
{
  CGAL::IO::set_ascii_mode(is);
  int n;
  is >> n;
  std::cout << "Reading " << n << " points" << std::endl;
  typename Weighted_point::Point p;
  for( ; n>0 ; n--)
    {
      if(is >> p) {
        L.push_back(Weighted_point(p,5*(n/10)));
      }
    }
  std::cout << "Points inserted" << std::endl;
  return true;
}

template <class AS>
void
_test_weighted_alpha_shape_3()
{
  typedef AS                                             Alpha_shape_3;

  typedef typename AS::Triangulation::Weighted_point     Weighted_point;
  typedef typename AS::Triangulation::Bare_point         Bare_point;

  typedef typename AS::Alpha_iterator                    Alpha_iterator;

  std::list<Weighted_point> L;
  bool verbose = false;

  // first a known small case
  // four groups of sphere :
  // - two groups of 4 intersecting spheres
  // - one group of three intersecting sphere
  // - one group of two intersecting sphere
  // the four groups are disjoint
  // Check specially the  $0$-shape which is the nerve of the union
  L.push_back(Weighted_point(Bare_point(0.,0.,0.), 4));
  L.push_back(Weighted_point(Bare_point(2.,2.,0.), 4));
  L.push_back(Weighted_point(Bare_point(2.,0.,2.), 4));
  L.push_back(Weighted_point(Bare_point(0.,2.,2.), 4));

  L.push_back(Weighted_point(Bare_point(10.,0.,0.), 4));
  L.push_back(Weighted_point(Bare_point(12.,2.,0.), 4));
  L.push_back(Weighted_point(Bare_point(12.,0.,2.), 4));
  L.push_back(Weighted_point(Bare_point(10.,2.,2.), 4));

  L.push_back(Weighted_point(Bare_point(20.,0.,0.), 3));
  L.push_back(Weighted_point(Bare_point(22.,0.,0.), 3));
  L.push_back(Weighted_point(Bare_point(20.,2.,0.), 3));

  L.push_back(Weighted_point(Bare_point(30.,0.,0.), 3));
  L.push_back(Weighted_point(Bare_point(32.,0.,0.), 3));


  Alpha_shape_3 A( L.begin(), L.end(), 0, Alpha_shape_3::REGULARIZED);
  if(verbose) show_alpha_values(A);
  A.set_alpha(0.);
  count_faces(A, verbose);
  test_filtration(A,verbose);

  A.set_mode(Alpha_shape_3::GENERAL);
  if(verbose) show_alpha_values(A);
  A.set_alpha(0.);
  count_faces(A, verbose);
  test_filtration(A,verbose);

  A.set_mode(Alpha_shape_3::REGULARIZED);
  assert(A.number_of_solid_components(0.) == 2);
  assert(A.number_of_solid_components(*(A.find_optimal_alpha(2))) <= 2);
  assert(A.number_of_solid_components(*(A.find_optimal_alpha(1))) == 1);
  if(verbose) show_alpha_values(A);
  A.set_alpha(0.);
  count_faces(A, verbose);
  test_filtration(A,verbose);

  // I add a test for CH4
  //  This test detected a bug in make_alpha_shape() cominig from clear()
  if(verbose) std::cerr << "test for CH4" << std::endl;
  L.clear();
  L.push_back(Weighted_point(Bare_point(-1.,-1.,-1.), 1.));
  L.push_back(Weighted_point(Bare_point(-1., 1., 1.), 1.));
  L.push_back(Weighted_point(Bare_point( 1.,-1., 1.), 1.));
  L.push_back(Weighted_point(Bare_point( 1., 1.,-1.), 1.));
  L.push_back(Weighted_point(Bare_point( 0., 0., 0.), 3.));
  Alpha_shape_3 A2( L.begin(), L.end(), 0, Alpha_shape_3::GENERAL);
  A2.set_alpha(0.);
  count_faces(A2, verbose);
  if(verbose) show_triangulation(A2);
  if(verbose) A2.print_alphas();
  if(verbose) A2.print_maps();
  assert(A2.get_mode()==Alpha_shape_3::GENERAL);
  if(verbose) std::cout << std::endl;
  if(verbose) std::cout << "test CH4 through make_alpha_shape" << std::endl;
  Alpha_shape_3 a2(0, Alpha_shape_3::GENERAL);
  a2.make_alpha_shape( L.begin(), L.end());
  count_faces(a2, verbose);
  if (verbose) show_triangulation(a2);
  if(verbose) a2.print_alphas();
  if(verbose) a2.print_maps();
  assert(a2.get_mode()==Alpha_shape_3::GENERAL);
  if(verbose) std::cout << std::endl;

  // test a bigger Alpha_shapes
  A.clear();
  L.clear();
  std::ifstream is("./data/fin", std::ios::in);
  assert(is);
  file_input(is,L);
  A.set_mode(Alpha_shape_3::REGULARIZED);
  std::size_t  n = A.make_alpha_shape(L.begin(), L.end());
  std::cout << "Alpha Shape computed :" << n  << " points" << std::endl;
  std::cout << " test number_of_components - find_optimal_alpha "
            <<   std::endl;
  A.set_alpha(*A.find_optimal_alpha(2));
  assert( A.number_of_solid_components() <= 2);

  Alpha_iterator opt = A.find_optimal_alpha(1);
  Alpha_iterator previous = opt; --previous;
  if(verbose) {
    std::cerr << " alpha optimal for 1 component  = " << *opt
              << "nb of components " << A.number_of_solid_components(*opt)
              << std::endl;
    std::cerr << " previous        " << *previous
              << "nb of components "
              << A.number_of_solid_components(*previous) << std::endl;
    std::cerr << "alpha_solid " << A.find_alpha_solid() << std::endl;
  }
  assert (A.number_of_solid_components(*opt) == 1);
  assert (A.number_of_solid_components(*previous) > 1
          || *opt ==  A.find_alpha_solid());
}



#endif //CGAL_TEST_WEIGHTED_ALPHA_SHAPE_3_H
