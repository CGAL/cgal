// Copyright (c) 1998-2003,2009,2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_TEST_CLS_ALPHA_SHAPE_3_H
#define CGAL_TEST_CLS_ALPHA_SHAPE_3_H

#include <list>
#include <fstream>
#include "_count_alpha_periodic_3.h"

template <class Point>
bool
file_input(std::ifstream& is, std::list<Point>& L, int nb=0)
{
  CGAL::IO::set_ascii_mode(is);
  int n;
  is >> n;
  if(nb != 0 && nb <= n) n=nb;
  std::cout << "Reading " << n << " points" << std::endl;
  Point p;
  for(; n>0; n--) {
    is >> p;
    L.push_back(p);
  }
  std::cout << "Points inserted" << std::endl;
  return true;
}

template <class AS>
void
_test_cls_alpha_shape_3()
{
  typedef AS                          Alpha_shape_3;
  typedef typename AS::Triangulation  Triangulation;

  typedef typename AS::FT FT;
  typedef typename AS::Point Point;
  typedef typename AS::Alpha_iterator Alpha_iterator;
  typedef typename Triangulation::Iso_cuboid Iso_cuboid;

  std::list<Point> L;
  bool verbose = false;

  // test a bigger alpha_shapes
  std::ifstream is("data/P3DT3_alpha_shape_test.pts", std::ios::in);
  assert(is);
  file_input(is, L);

  Iso_cuboid domain(FT(-0.1),FT(-0.1),FT(-0.1), FT(0.2),FT(0.2),FT(0.2));
  Triangulation dt(domain);
  std::cout << "Create triangulation with domain " << domain << std::endl;
  std::size_t n = dt.insert(L.begin(),L.end());

  Alpha_shape_3 a1(dt,0, Alpha_shape_3::REGULARIZED);
  if(verbose)
    show_alpha_values(a1);

  std::cout << "Alpha Shape computed: " << n << " points" << std::endl;
  std::cout << " test number_of_components - find_optimal_alpha "<< std::endl;

  Alpha_iterator opt = a1.find_optimal_alpha(1);
  Alpha_iterator previous = opt;
  --previous;

  if(verbose) {
    std::cerr << " optimal de 1 " << *opt
              << "nb of components " << a1.number_of_solid_components(*opt)
              << std::endl;
    std::cerr << " previous " << *previous
              << "nb of components "
              << a1.number_of_solid_components(*previous) << std::endl;
  }
  assert(a1.number_of_solid_components(*opt) == 1);
  assert(a1.number_of_solid_components(*previous) > 1);
}

template <class AS>
void
_test_cls_alpha_shape_3_exact()
{
  typedef AS                          Alpha_shape_3;
  typedef typename AS::Triangulation  Triangulation;

  typedef typename AS::FT FT;
  typedef typename AS::Point Point;
  typedef typename AS::Alpha_iterator Alpha_iterator;

  bool verbose = false;
  std::vector<Point> L, Lc;

  L.push_back(Point(0, 0, 0));
  L.push_back(Point(1, 2, 3));
  Lc.push_back(Point(0, 0, 0));

  FT a(1), b(2), c(3), d(8);
  Lc.push_back(Point(a/d, b/d, c/d));

  Triangulation dt1(L.begin(), L.end(),
                    typename Triangulation::Iso_cuboid(0,0,0, 8,8,8));
  Alpha_shape_3 a1(dt1, 0, Alpha_shape_3::REGULARIZED);

  if(verbose)
    show_triangulation(a1);
  assert(a1.number_of_alphas() == 6);

  std::cout << "REGULARIZED mode" << std::endl;
  std::cout << "test_classify_and_iterators" << std::endl;

  Alpha_iterator alpha_it = a1.alpha_begin();
  for(;alpha_it != a1.alpha_end();alpha_it++) {
    a1.set_alpha(*alpha_it);
    if(verbose) {
      std::cerr << std::endl;
      std::cerr << "alpha value " << * alpha_it << std::endl;
      std::cerr << std::endl;
    }
    count_faces(a1, verbose);
  }

  if(verbose) show_alpha_values(a1);
  if(verbose) a1.print_maps();
  if(verbose) a1.print_alphas();

  std::cout << "test filtration " << std::endl;
  test_filtration(a1, verbose);

  a1.set_mode(Alpha_shape_3::GENERAL);
  std::cout << "GENERAL mode" << std::endl;
  if(verbose) std::cerr << "GENERAL mode" << std::endl;
  if(verbose) show_alpha_values(a1);
  if(verbose) a1.print_maps();
  if(verbose) a1.print_alphas();
  assert(a1.number_of_alphas() == 21);
  std::cout << "test_classify_and_iterators" << std::endl;

  for(alpha_it = a1.alpha_begin();alpha_it!=a1.alpha_end();alpha_it++) {
    a1.set_alpha(*alpha_it);
    if(verbose) {
      std::cerr << std::endl;
      std::cerr << "alpha value " << * alpha_it << std::endl;
    }
    count_faces(a1, verbose);
  }
  std::cout << "test filtration " << std::endl;
  test_filtration(a1, verbose);

  // alpha values 0
  //              21
  //              25.4444
  //              26
  //              29
  //              30.4444
  //              30.7959

  a1.set_mode(Alpha_shape_3::REGULARIZED);
  std::cout << "REGULARIZED mode" << std::endl;
  std::cout << "test number_of_components - find_optimal_alpha " << std::endl;

  if(verbose) {
    std::cerr << std::endl;
    std::cerr << "REGULARIZED mode" << std::endl;
    a1.print_maps();
    a1.print_alphas();

    for(alpha_it = a1.alpha_begin();alpha_it!=a1.alpha_end();alpha_it++)
      std::cerr << "alpha " << *alpha_it << "\t"
                << "number of solid componenst "
                << a1.number_of_solid_components(*alpha_it) << std::endl;
  }

  // alpha  21       number of solid components 54
  // alpha  25.4444  number of solid components 54
  // alpha  26       number of solid components 27
  // alpha  29       number of solid components  9
  // alpha  30.4444  number of solid components  3
  // alpha  30.7949  number of solid components  1

  assert(*(a1.find_optimal_alpha(1)) == a1.get_nth_alpha(6));
  assert(*(a1.find_optimal_alpha(5)) == a1.get_nth_alpha(5));
  assert(a1.number_of_solid_components(*(a1.find_optimal_alpha(1))) == 1);
  assert(a1.number_of_solid_components(*(a1.find_optimal_alpha(4))) == 3);
  assert(a1.get_nth_alpha(3) == *(a1.alpha_lower_bound(26)));
  assert(a1.get_nth_alpha(4) == *(a1.alpha_upper_bound(26)));
  test_filtration(a1, verbose);

  std::cout << std::endl;
  std::cout << "test additional creators and set mode" << std::endl;

  Triangulation dt2(Lc.begin(), Lc.end());
  Alpha_shape_3 a2(dt2, 0, Alpha_shape_3::REGULARIZED);
  if(verbose) show_alpha_values(a2);
  if(verbose) a2.print_maps();
  if(verbose) a2.print_alphas();
  assert(a2.number_of_alphas() == 6);

  Triangulation dt3(Lc.begin(), Lc.end());
  Alpha_shape_3 a3(dt3, 0, Alpha_shape_3::GENERAL);
  assert(a3.number_of_alphas() == 21);

  Triangulation dt4(Lc.begin(), Lc.end());
  Alpha_shape_3 a4(dt4);
  assert(a4.number_of_alphas() == 6);

  Alpha_shape_3 a5(0, Alpha_shape_3::GENERAL);
  a5.make_alpha_shape(Lc.begin(), Lc.end());
  assert(a5.number_of_alphas() == 21);

  a1.set_mode(Alpha_shape_3::REGULARIZED);
  assert(a1.number_of_alphas() == 6);

  a1.set_mode(Alpha_shape_3::GENERAL);
  assert(a1.number_of_alphas() == 21);
}

#endif // CGAL_TEST_CLS_ALPHA_SHAPE_3_H
