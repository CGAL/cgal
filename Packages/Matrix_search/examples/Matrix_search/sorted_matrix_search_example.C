// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : sorted_matrix_search_example.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : fjsearch.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Sorted matrix search: Example Program
// ============================================================================

#ifndef CGAL_RANDOM_H
#include <CGAL/Random.h>
#endif // CGAL_RANDOM_H
#ifndef CGAL_FUNCTION_OBJECTS_H
#include <CGAL/function_objects.h>
#endif // CGAL_FUNCTION_OBJECTS_H
#ifndef CGAL_CARTESIAN_MATRIX_H
#include <CGAL/Cartesian_matrix.h>
#endif // CGAL_CARTESIAN_MATRIX_H
#ifndef CGAL_SORTED_MATRIX_SEARCH_H
#include <CGAL/sorted_matrix_search.h>
#endif // CGAL_SORTED_MATRIX_SEARCH_H
#ifndef CGAL_PROTECT_VECTOR_H
#include <vector.h>
#define CGAL_PROTECT_VECTOR_H
#endif // CGAL_PROTECT_VECTOR_H

int main() {

  typedef int                              Value;
  typedef vector< Value >                  Vector;
  typedef Vector::iterator                 Value_iterator;
  typedef vector< Vector >                 Vector_cont;
  typedef CGAL_Cartesian_matrix<
    plus< int >,
    Value_iterator,
    Value_iterator >                       Matrix;

  // set of vectors the matrices are build from:
  Vector_cont vectors;

  // generate a random vector and sort it:
  Vector a;
  int i;
  cout << "a = ( ";
  for ( i = 0; i < 5; ++i) {
    a.push_back( CGAL_random( 100));
    cout << a.back() << " ";
  }
  cout << ")" << endl;
  sort( a.begin(), a.end(), less< Value >());

  // build a cartesian from a:
  Matrix M( a.begin(), a.end(), a.begin(), a.end());

  // search an upper bound for max(a):
  Value bound( a[4]);
  Value upper_bound(
    CGAL_sorted_matrix_search(
      &M,
      &M + 1,
      CGAL_sorted_matrix_search_traits_adaptor(
        bind2nd( greater_equal< Value >(), bound),
        M)));
  cout << "upper bound for " << bound << " is "
       << upper_bound << endl;

} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

