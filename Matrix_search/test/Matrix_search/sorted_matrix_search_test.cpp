// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#include <CGAL/Random.h>
#include <CGAL/Cartesian_matrix.h>
#include <CGAL/sorted_matrix_search.h>
#include <vector>
#include <cstdlib>

template < class Matrix_iterator, class Value >
Value
compute_upper_bound( Matrix_iterator f,
                     Matrix_iterator l,
                     Value b,
                     Value max_bound)
// computes the smallest element of one of the matrices
// from the range [f,l[ that is greater than b and smaller
// or equal to max_bound
// brute force in O( sum_{m \in [f,l[} m.x_dimension() * m.y_dimension() )
//
{
  Value best( max_bound);
  Matrix_iterator k( f);
  while ( k != l) {
    for ( int i = 0; i < (*k).number_of_columns(); ++i)
      for ( int j = 0; j < (*k).number_of_rows(); ++j)
        if ( (*k)( i, j) >= b && (*k)( i, j) < best)
          best = (*k)( i, j);
    ++k;
  }
  return best;
} // compute_upper_bound( f, l, b, max)
using std::vector;
using std::plus;
using std::sort;
using std::less;
using std::greater_equal;
using std::max;
using std::cerr;
using std::cout;
using std::endl;
using std::atoi;
using std::exit;
using CGAL::Cartesian_matrix;
using CGAL::Random;
using CGAL::get_default_random;
using CGAL::sorted_matrix_search;
using CGAL::sorted_matrix_search_traits_adaptor;

typedef int                              Value;
typedef vector< Value >                  Vector;
typedef Vector::iterator                 Value_iterator;
typedef vector< Vector >                 Vector_cont;
typedef Vector_cont::iterator            Vector_iterator;
typedef Cartesian_matrix<
  plus< int >,
  Value_iterator,
  Value_iterator >                       Matrix;
typedef vector< Matrix >                 Matrix_cont;



int
main( int argc, char* argv[])
{
  // seed for random number generator:
  int random_seed;
  // number of matrices:
  int num;
  // dimension of matrices to generate:
  int dim;
  // bound for matrix entries:
  int upper_entry_bound;
  // set of matrices to search:
  Matrix_cont matrices;
  // set of vectors the matrices are build from:
  Vector_cont vectors;
  // handle command line arguments:
  if ( argc < 4 ||
       (num = atoi(argv[1])) <= 0 ||
       (dim = atoi(argv[2])) <= 0 ||
       (upper_entry_bound = atoi(argv[3])) <= 0) {
    cerr << "usage: " << argv[0] <<
      " num dim upper_entry_bound [random_seed]" << endl;
    exit(1);
  }
  if ( argc < 5 || (random_seed = atoi(argv[4])) <= 0) {
  
    #ifdef OUTPUT
    cerr << "No random seed specified - generating it" << endl;
    #endif
  
    // generate random seed
    random_seed = get_default_random().get_int( 0, (1 << 30));
  }
  else
    random_seed = atoi(argv[4]);
  
  // define random source:
  Random r( random_seed);
  
  #ifdef OUTPUT
  cout << "random seed is " << random_seed << endl;
  #endif
  // maximum entry of all matrices:
  Value max_entry( -1);
  
  for ( int k = 0; k < num; ++k) {
    // generate two vectors a and b to build the cartesian
    // matrix from
    Vector a, b;
    assert( a.size() == 0 && b.size() == 0);
  
    // fill a and b with random values and sort them:
    for ( int i = 0; i < dim; ++i) {
      a.push_back( r( upper_entry_bound));
      b.push_back( r( upper_entry_bound));
    }
    // to test some non-quadratic matrices:
    // for ( i = 0; i < dim / 5; ++i) {
    //   b.push_back( r());
    // }
  
    sort( a.begin(), a.end(), less< Value >());
    sort( b.begin(), b.end(), less< Value >());
  
    /*
    cout << "a = (";
    for ( Vector_iterator pp( a.begin());
    pp != a.end();
    ++pp)
      cout << (*pp) << ", ";
    cout << ")\nb = (";
    for ( Vector_iterator pq( a.begin());
    pq != a.end();
    ++pq)
      cout << (*pq) << ", ";
    cout << ")" << endl;
    */
  
    // evt. update max_entry:
    max_entry = max( a[dim - 1] + b[dim - 1], max_entry);
  
    // keep both vectors:
    vectors.push_back( a);
    vectors.push_back( b);
  } // for ( int k = 0; k < num; ++k)
  
  
  // construct matrices:
  for ( Vector_iterator i( vectors.begin());
        i != vectors.end();
        i += 2)
    {
      Vector_iterator j = i + 1;
      matrices.push_back(
        Matrix( (*i).begin(), (*i).end(),
                (*j).begin(), (*j).end()));
    }
  // search lower bound for a random value v in matrices
  Value bound;
  // assure there is any feasible value in m:
  do
    bound = r.get_int( 0, 2 * upper_entry_bound);
  while ( bound > max_entry);
  
  #ifdef OUTPUT
  cout << "searching upper bound for " << bound << endl;
  #endif
  Value u(
    sorted_matrix_search(
      matrices.begin(),
      matrices.end(),
      sorted_matrix_search_traits_adaptor(
        boost::bind( greater_equal< Value >(), _1, bound),
        *(matrices.begin()))));
  
  #ifdef OUTPUT
  cout << "************* DONE *************\nresult: "
       << u << "\n********************************" << endl;
  CGAL_optimisation_assertion(
    u == compute_upper_bound(
      matrices.begin(), matrices.end(), bound, max_entry));
  #else
  Value brute_force(
    compute_upper_bound(
      matrices.begin(), matrices.end(), bound, max_entry));
  if ( u != brute_force)
    cerr << "\nerror: num( " << num << "), dim ( "
         << dim << "), upper_entry_bound( " << upper_entry_bound
         << ")\nrandom_seed( " << random_seed << ")"
         << "\nresult was " << u << "\ntrivial algorithm gives "
         << brute_force << endl;
  #endif

  return 0;
}
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

