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
// file          : monotone_matrix_search.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Monotone Matrix Search
// ============================================================================

#if ! (CGAL_MONOTONE_MATRIX_SEARCH_H)
#define CGAL_MONOTONE_MATRIX_SEARCH_H 1

#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#include <CGAL/optimisation_assertions.h>
#endif // CGAL_OPTIMISATION_ASSERTIONS_H
#ifndef CGAL_PROTECT_VECTOR_H
#include <vector.h>
#define CGAL_PROTECT_VECTOR_H
#endif // CGAL_PROTECT_VECTOR_H
#ifndef CGAL_PROTECT_FUNCTION_H
#include <function.h>
#define CGAL_PROTECT_FUNCTION_H
#endif // CGAL_PROTECT_FUNCTION_H

template < class Matrix, class RandomAccessIterator >
inline
void
CGAL_monotone_matrix_search(
  const Matrix& M,
  RandomAccessIterator t)
{
  typedef typename Matrix::Value V;
  CGAL_monotone_matrix_search( M, t, less< V >());
} // CGAL_monotone_matrix_search( M, t)
template < class Matrix,
           class RandomAccessIterator,
           class Compare_strictly >
void
CGAL_monotone_matrix_search(
  const Matrix& M,
  RandomAccessIterator t,
  const Compare_strictly& compare_strictly)
// Matrix has to define:
//  o operator()( int, int) [access]
//  o int number_of_columns(), int number_of_rows()
//  o replace_column( int o, int n)
//  o shrink_to_quadratic_size()
//  o Matrix* extract_all_even_rows()
//
// Precondition: M is totally monotone
//  M.number_of_rows() > 1 and
// RandomAccessIterator has value type int
//
// writes to t the positions (columns)
// of the row maxima of M
{
  // divide
  // ------
  // get even rows of M:
  #ifdef CGAL_MON_SEARCH_TRACE
  cerr << "construct new matrix" << endl;
  #endif
  Matrix* M_new( M.extract_all_even_rows());
  CGAL_optimisation_assertion(
    M_new->number_of_columns() == M.number_of_columns());
  CGAL_optimisation_assertion(
    M_new->number_of_rows() == 0 ||
      M_new->number_of_rows() == ( M.number_of_rows() + 1) >> 1);
  
  #ifdef CGAL_MON_SEARCH_TRACE
  {
    for ( int i1( 0); i1 < M_new->number_of_rows(); ++i1) {
      for ( int i2( 0); i2 < M_new->number_of_columns(); ++i2) {
        cerr << i1 << " - " << i2 << endl;
        cout.width( 4);
        cout << (*M_new)( i1, i2) << "  ";
      }
      cout << endl;
    }
    cout << "----------------------" << endl;
  }
  #endif

  // reduce M_new to a quadratic matrix:
  #ifdef CGAL_MON_SEARCH_TRACE
  cerr << "reduce" << endl;
  #endif
  
  // table to store the reduction permutation:
  // (incl. sentinel)
  int* reduction_table( new int[ M_new->number_of_rows() + 1]);
  
  if ( M_new->number_of_rows() < M_new->number_of_columns()) {
    // set sentinel:
    reduction_table[M_new->number_of_rows()] =
      M.number_of_columns() - 1;
    CGAL__reduce_matrix( *M_new, reduction_table, compare_strictly);
    CGAL_optimisation_assertion(
      M_new->number_of_columns() == M_new->number_of_rows());
  
    #ifdef CGAL_MON_SEARCH_TRACE
    {
      int i1;
      for ( i1 = 0; i1 < M_new->number_of_rows(); ++i1) {
        for ( int i2( 0); i2 < M_new->number_of_columns(); ++i2) {
          cout.width( 4);
          cout << (*M_new)( i1, i2) << "  ";
        }
        cout << endl;
      }
      cout << "----------------------\n reduction table:" << endl;
      for ( i1 = 0; i1 < M_new->number_of_rows(); ++i1) {
        cout << reduction_table[i1] << ", ";
      }
      cout << "\n----------------------" << endl;
    }
    #endif
  } // if ( M_new->number_of_rows() < M_new->number_of_columns())
  else {
    // no reduction -> reduction_table is identity table:
    for ( int i1( 0); i1 < M_new->number_of_columns(); ++i1)
      reduction_table[i1] = i1;
    // set sentinel:
    reduction_table[M_new->number_of_columns()] =
      M_new->number_of_columns() - 1;
  }
  
  

  // recursion:
  
  CGAL_optimisation_assertion(
    M_new->number_of_rows() >= M_new->number_of_columns());
  
  // table to store the rmax values of M_new:
  // (incl. sentinel)
  int* t_new( new int[M_new->number_of_rows() + 1]);
  t_new[M_new->number_of_rows()] = M_new->number_of_columns();
  
  #ifdef CGAL_MON_SEARCH_TRACE
  cerr << "recursive call" << endl;
  #endif
  if ( M_new->number_of_rows() == 1)
    // recursion anchor:
    // we have just one element ==> no choice
    t_new[0] = 0;
  else
    CGAL_monotone_matrix_search( *M_new, t_new);
  
  #ifdef CGAL_MON_SEARCH_TRACE
  cerr << "maximum entries:\n";
  int i = 0;
  while ( i < M_new->number_of_rows())
    cerr << t_new[i++] << "  ";
  cerr << endl;
  #endif

  // and conquer
  // -----------
  #ifdef CGAL_MON_SEARCH_TRACE
  {
    cerr << "find maxima in odd rows" << endl;
    for ( int i1( 0); i1 < M.number_of_rows(); ++i1) {
      for ( int i2( 0); i2 < M.number_of_columns(); ++i2) {
        cout.width( 4);
        cout << M( i1, i2) << "  ";
      }
      cout << endl;
    }
    cout << "----------------------" << endl;
  }
  #endif
  
  int j( 0);       // actual index in t
  int j_new( 0);   // actual index in t_new
  do {
    // even row ==> we know
    *(t+j) = reduction_table[t_new[j_new++]];
    #ifdef CGAL_MON_SEARCH_TRACE
    cerr << " # maximum of row " << j << " was at " << *(t+j) << endl;
    #endif
    if ( ++j >= M.number_of_rows())
      break;
  
    // odd row
    // search *(t+j) between *(t+j-1) and t_new[j_new]:
    #ifdef CGAL_MON_SEARCH_TRACE
    cerr << "search row " << j << " between " <<
      *(t+j-1) << " and " <<
      reduction_table[t_new[j_new]] << endl;
    #endif
    *(t+j) = reduction_table[t_new[j_new]];
    int j_tmp( *(t+j-1));
    while ( j_tmp < reduction_table[t_new[j_new]]) {
      if ( compare_strictly( M( j, t[j]), M( j, j_tmp)))
        *(t+j) = j_tmp++;
      else
        ++j_tmp;
    }
  } while ( ++j < M.number_of_rows());

  delete M_new;
  delete[] t_new;
  delete[] reduction_table;

} // CGAL_monotone_matrix_search( M, t)


template < class Matrix,
           class RandomAccessIterator,
           class Compare_strictly >
void
CGAL__reduce_matrix(
  Matrix& M,
  RandomAccessIterator t,
  const Compare_strictly& compare_strictly)
// Matrix has to define:
//  o operator()( int, int) [access]
//  o int number_of_columns(), int number_of_rows()
//  o replace_column( int o, int n)
//  o shrink_to_quadratic_size()
//  o Matrix* extract_all_even_rows()
//
// Precondition: M is totally monotone
// reduces M, i.e. deletes some columns that
// do not contain the maximum value of any row
// such that M becomes quadratic
// and returns for each column of the resulting
// matrix its column index in the original matrix
{
  CGAL_optimisation_precondition(
    M.number_of_columns() >= M.number_of_rows());
  // active columns are 0, ..., j1, j2, ..., M.x_dim()-1
  int j1( 0), j2( 1);
  *t = 0;
  while ( j2 - j1 < M.number_of_columns() - M.number_of_rows() + 1) {
    if ( compare_strictly( M( j1, j1), M( j1, j2))) {
      // delete column j1
      if ( j1 > 0)
        --j1;
      else {
        M.replace_column( 0, j2);
        *t = j2++;
      }
    }
    else {
      if ( j1 < M.number_of_rows() - 1) {
        // proceed
        M.replace_column( ++j1, j2);
        *(t+j1) = j2;
      }
      // else delete column j2
      ++j2;
    }
  } // while ( j2 - j1 <
    //         M.number_of_columns() - M.number_of_rows() + 1)
  
  // M.number_of_columns() - M.number_of_rows() columns
  // have been deleted, now move columns
  // j2 .. M.number_of_columns()-1 to the first part
  while ( j1 < M.number_of_rows() - 1) {
    CGAL_optimisation_assertion( j2 < M.number_of_columns());
    M.replace_column( ++j1, j2);
    *(t+j1) = j2++;
  }
  
  M.shrink_to_quadratic_size();
} // CGAL__reduce_matrix( M, t)

#endif // ! (CGAL_MONOTONE_MATRIX_SEARCH_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

