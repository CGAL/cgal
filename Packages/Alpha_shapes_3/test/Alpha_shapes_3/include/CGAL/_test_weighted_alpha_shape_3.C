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

#ifndef CGAL_TEST_WEIGHTED_ALPHA_SHAPE_3_C
#define CGAL_TEST_WEIGHTED_ALPHA_SHAPE_3_C

#include <list>
#include <fstream>
#include "_count_alpha.C"

template <class Weighted_point>
bool
file_input(std::ifstream& is, std::list<Weighted_point>& L)
{
  CGAL::set_ascii_mode(is);
  int n;
  is >> n;
  std::cout << "Reading " << n << " points" << std::endl;
  typename Weighted_point::Point p;
  for( ; n>0 ; n--)
    {
      is >> p;
      L.push_back(Weighted_point(p,5*(n/10)));
    }
  std::cout << "Points inserted" << std::endl;
  return true;
}

template <class AS>
void
_test_weighted_alpha_shape_3()
{
  typedef AS                 Alpha_shape_3;

  //typedef typename AS::Point Point;
  typedef typename AS::Gt::Weighted_point Weighted_point;
  typedef typename AS::Gt::Bare_point     Bare_point;
  typedef typename AS::Alpha_iterator Alpha_iterator;
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

  A.set_mode(Alpha_shape_3::GENERAL);
  if(verbose) show_alpha_values(A);
  A.set_alpha(0.);
  count_faces(A, verbose);

  assert(A.number_of_solid_components(0.) == 2);
  assert(A.number_of_solid_components(*(A.find_optimal_alpha(2))) <= 2);
  assert(A.number_of_solid_components(*(A.find_optimal_alpha(1))) == 1);

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
  A.set_mode(Alpha_shape_3::GENERAL);
  int  n = A.make_alpha_shape(L.begin(), L.end());
  std::cout << "Alpha Shape computed :" << n  << " points" << std::endl;
  std::cout << " test number_of_components - find_optimal_alpha "
	    <<   std::endl;
  A.set_alpha(*A.find_optimal_alpha(2));
  //show_alpha_values(A);
  // std::cerr << "optimal alpha " << *A.find_optimal_alpha(2) << std::endl;
  assert( A.number_of_solid_components() <= 2);
  assert(A.number_of_solid_components(*(A.find_optimal_alpha(1))) == 1);

}



#endif //CGAL_TEST_WEIGHTED_ALPHA_SHAPE_3_C
