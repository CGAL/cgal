// Copyright (c) 1998-2003  INRIA Sophia-Antipolis (France).
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
// $Source$ 
// $Revision$ 
// $Date$
// $Name$
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_TEST_CLS_ALPHA_SHAPE_3_C
#define CGAL_TEST_CLS_ALPHA_SHAPE_3_C

#include <list>
#include <fstream>
#include "_count_alpha.C"


template <class Point>
bool
file_input(std::ifstream& is, std::list<Point>& L)
{
  CGAL::set_ascii_mode(is);
  int n;
  is >> n;
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
_test_cls_alpha_shape_3(const AS &)
{
  typedef AS                 Alpha_shape_3;

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
  
  Alpha_shape_3 A( L.begin(), L.end(), 0, Alpha_shape_3::REGULARIZED);
  if(verbose) show_triangulation(A);
  assert(A.number_of_alphas() == 3 );
  std::cout << "test_classify_and_iterators in regularised mode" 
	    << std::endl;
  Alpha_iterator alpha_it = A.alpha_begin();
  for (;alpha_it != A.alpha_end();alpha_it++){
    A.set_alpha(*alpha_it);
    if (verbose) {
      std::cerr << std::endl;
      std::cerr << "alpha value " << * alpha_it << std::endl;
    }
    count_faces(A, verbose);
  }

  A.set_mode(Alpha_shape_3::GENERAL);
  if(verbose) show_alpha_values(A);
  if(verbose) A.print_maps();
  if(verbose) A.print_alphas();
  assert(A.number_of_alphas() == 8) ;
	 std::cout << "test_classify_and_iterators in general mode" 
	    << std::endl;
 
  for(alpha_it = A.alpha_begin();alpha_it!=A.alpha_end();alpha_it++){
    A.set_alpha(*alpha_it);
    if (verbose) {
      std::cerr << std::endl;
      std::cerr << "alpha value " << * alpha_it << std::endl;
    }
    count_faces(A, verbose);
  }

  // alpha values 1
  //              1.0625
  //              1.38942
  //              2
  //              2.00694
  //              2.66667
  //              3
  //              4.35417

  std::cout << "test number_of_components - find_optimal_alpha "
	    <<  std::endl;
  if (verbose) {
    std::cerr << std::endl;
    for(alpha_it = A.alpha_begin();alpha_it!=A.alpha_end();alpha_it++)
      std::cerr << "alpha  " << *alpha_it << "\t" 
		<< "number of solid componenst " 
		<< A.number_of_solid_components(*alpha_it) <<  std::endl;
  }
  assert( *(A.find_optimal_alpha(1)) == A.get_nth_alpha(8));
  assert( *(A.find_optimal_alpha(2)) == A.get_nth_alpha(6));
  assert (A.number_of_solid_components(*(A.find_optimal_alpha(2))) == 2);
  assert (A.number_of_solid_components(*(A.find_optimal_alpha(1))) == 1);
  assert (A.get_nth_alpha(1) == *(A.alpha_lower_bound(1)));
  assert (A.get_nth_alpha(5) == *(A.alpha_upper_bound(2)));

  // test a bigger alpha_shapes
  A.clear();
  L.clear();
  std::ifstream is("./data/fin", std::ios::in);
  assert(is);
  file_input(is,L);
  A.set_mode(Alpha_shape_3::GENERAL);
  int  n = A.make_alpha_shape(L.begin(), L.end());
  std::cout << "Alpha Shape computed :"  << n  << " points" << std::endl;
  A.set_alpha(*A.find_optimal_alpha(2));
  std::cout << " test number_of_components - find_optimal_alpha "<< std::endl;
  assert( A.number_of_solid_components() == 2);
}

#endif //  CGAL_TEST_CLS_ALPHA_SHAPE_3_C


