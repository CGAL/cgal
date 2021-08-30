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

#ifndef CGAL_TEST_CLS_ALPHA_SHAPE_3_H
#define CGAL_TEST_CLS_ALPHA_SHAPE_3_H

#include <list>
#include <fstream>
#include "_count_alpha.h"


template <class Point>
bool
file_input(std::ifstream& is, std::list<Point>& L, int nb=0)
{
  CGAL::IO::set_ascii_mode(is);
  int n;
  is >> n;
  if (nb != 0 && nb <= n) n=nb;
  std::cout << "Reading " << n << " points" << std::endl;
  Point p;
  for( ; n>0 ; n--)    {
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

  typedef typename AS::Point Point;
  typedef typename AS::Alpha_iterator Alpha_iterator;


  std::list<Point> L;
  bool verbose = false;

  // first a known small case
  // a cube with one corner less and two small pyramids
  // on back and front face
  L.push_back(Point(0.,0.,0.));
  L.push_back(Point(2.,0.,0.));
  L.push_back(Point(0.,2.,0.));
  L.push_back(Point(0.,0.,2.));
  L.push_back(Point(2.,0.,2.));
  L.push_back(Point(0.,2.,2.));
  L.push_back(Point(2.,2.,0));
  L.push_back(Point(1.,-1.5 ,1.));
  L.push_back(Point(1.,3.5 ,1.));


  Alpha_shape_3 a1( L.begin(), L.end(), 0, Alpha_shape_3::REGULARIZED);
  if(verbose) show_triangulation(a1);
  assert(a1.number_of_alphas() == 3 );

  std::cout << "REGULARIZED  mode" << std::endl;
  std::cout << "test_classify_and_iterators"
            << std::endl;
  Alpha_iterator alpha_it = a1.alpha_begin();
  for (;alpha_it != a1.alpha_end();alpha_it++){
    a1.set_alpha(*alpha_it);
    if (verbose) {
      std::cerr << std::endl;
      std::cerr << "alpha value " << * alpha_it << std::endl;
      std::cerr << std::endl;
     }
    count_faces(a1, verbose);
  }

  if(verbose) show_alpha_values(a1);
  if (verbose)     a1.print_maps();
  if (verbose)   a1.print_alphas();
  std::cout << "test filtration " << std::endl;
  test_filtration(a1,verbose);

  a1.set_mode(Alpha_shape_3::GENERAL);
  std::cout << "GENERAL mode" << std::endl;
  if(verbose) std::cerr << "GENERAL mode" << std::endl;
  if(verbose) show_alpha_values(a1);
  if(verbose) a1.print_maps();
  if(verbose) a1.print_alphas();
  assert(a1.number_of_alphas() == 9) ;
         std::cout << "test_classify_and_iterators"
            << std::endl;

  for(alpha_it = a1.alpha_begin();alpha_it!=a1.alpha_end();alpha_it++){
    a1.set_alpha(*alpha_it);
    if (verbose) {
      std::cerr << std::endl;
      CGAL::to_double(* alpha_it);
      std::cerr << "alpha value " << * alpha_it << std::endl;
    }
    count_faces(a1, verbose);
  }
  std::cout << "test filtration " << std::endl;
  test_filtration(a1,verbose);

  // alpha values 0
  //              1
  //              1.0625
  //              1.38942
  //              2
  //              2.00694
  //              2.66667
  //              3
  //              4.35417


  a1.set_mode(Alpha_shape_3::REGULARIZED);
  std::cout << "REGULARIZED mode" << std::endl;
  std::cout << "test number_of_components - find_optimal_alpha "
            <<  std::endl;
  if (verbose) {
    std::cerr << std::endl;
    std::cerr << "REGULARIZED mode" << std::endl;
    a1.print_maps();
    a1.print_alphas();

    for(alpha_it = a1.alpha_begin();alpha_it!=a1.alpha_end();alpha_it++)
      std::cerr << "alpha  " << *alpha_it << "\t"
                << "number of solid componenst "
                << a1.number_of_solid_components(*alpha_it) <<  std::endl;
  }

  // alpha  2.00694  number of solid components 2
  // alpha  3        number of solid component 1
  // alpha  4.35417  number of solid components 1


  assert( *(a1.find_optimal_alpha(1)) == a1.get_nth_alpha(2));
  assert( *(a1.find_optimal_alpha(2)) == a1.get_nth_alpha(1));
  assert (a1.number_of_solid_components(*(a1.find_optimal_alpha(2))) == 2);
  assert (a1.number_of_solid_components(*(a1.find_optimal_alpha(1))) == 1);
  assert (a1.get_nth_alpha(2) == *(a1.alpha_lower_bound(3)));
  assert (a1.get_nth_alpha(3) == *(a1.alpha_upper_bound(4)));
  test_filtration(a1,verbose);

  std::cout << std::endl;
  std::cout << "test additionnal creators and set mode" << std::endl;
  Triangulation dt2( L.begin(), L.end());
  Alpha_shape_3 a2(dt2, 0,  Alpha_shape_3::REGULARIZED);
  if(verbose) show_alpha_values(a2);
  if(verbose) a2.print_maps();
  if(verbose) a2.print_alphas();
  assert(a2.number_of_alphas() == 3) ;

  Triangulation dt3( L.begin(), L.end());
  Alpha_shape_3 a3(dt3, 0,  Alpha_shape_3::GENERAL);
  assert(a3.number_of_alphas() == 9) ;

  Triangulation dt4( L.begin(), L.end());
  Alpha_shape_3 a4(dt4);
  assert(a4.number_of_alphas() == 3) ;

  Alpha_shape_3 a5(0,  Alpha_shape_3::GENERAL);
  a5.make_alpha_shape(L.begin(), L.end());
  assert(a5.number_of_alphas() == 9) ;

  a1.set_mode(Alpha_shape_3::REGULARIZED);
  assert(a1.number_of_alphas() == 3) ;

  a1.set_mode(Alpha_shape_3::GENERAL);
  assert(a1.number_of_alphas() == 9) ;

// test a bigger alpha_shapes
  a1.clear();
  L.clear();
  std::ifstream is("./data/fin", std::ios::in);
  assert(is);
  file_input(is,L);
  a1.set_mode(Alpha_shape_3::REGULARIZED);
  std::size_t  n = a1.make_alpha_shape(L.begin(), L.end());
  if(verbose) show_alpha_values(a1);
  std::cout << "Alpha Shape computed :"  << n  << " points" << std::endl;
  std::cout << " test number_of_components - find_optimal_alpha "<< std::endl;
  Alpha_iterator opt = a1.find_optimal_alpha(1);
  Alpha_iterator previous = opt; --previous;
  if(verbose) {
    std::cerr << " optimal  de 1 " << *opt
              << "nb of componants " << a1.number_of_solid_components(*opt)
              << std::endl;
    std::cerr << " previous        " << *previous
              << "nb of componants "
              << a1.number_of_solid_components(*previous) << std::endl;
  }
  assert (a1.number_of_solid_components(*opt) == 1);
  assert (a1.number_of_solid_components(*previous) > 1);
}



#endif //  CGAL_TEST_CLS_ALPHA_SHAPE_3_H
