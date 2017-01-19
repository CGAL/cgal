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

#ifndef CGAL_MONOTONE_MATRIX_SEARCH_H
#define CGAL_MONOTONE_MATRIX_SEARCH_H 1

#include <CGAL/license/Matrix_search.h>


#include <CGAL/Optimisation/assertions.h>
#include <vector>
#include <functional>

namespace CGAL {
template < class Matrix, class RandomAccessIterator >
inline
void
monotone_matrix_search(
  const Matrix& M,
  RandomAccessIterator t)
{
  typedef typename Matrix::Value V;
  monotone_matrix_search( M, t, std::less< V >());
} // monotone_matrix_search( M, t)
template < class Matrix,
           class RandomAccessIterator,
           class Compare_strictly >
void
monotone_matrix_search(
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
  Matrix* M_new = M.extract_all_even_rows();
  CGAL_optimisation_assertion(
    M_new->number_of_columns() == M.number_of_columns());
  CGAL_optimisation_assertion(
    M_new->number_of_rows() == 0 ||
      M_new->number_of_rows() == ( M.number_of_rows() + 1) >> 1);
  

  // reduce M_new to a quadratic matrix:
  
  // table to store the reduction permutation:
  // (incl. sentinel)
  int* reduction_table = new int[ M_new->number_of_rows() + 1];
  
  if ( M_new->number_of_rows() < M_new->number_of_columns()) {
    // set sentinel:
    reduction_table[M_new->number_of_rows()] =
      M.number_of_columns() - 1;
    _reduce_matrix( *M_new, reduction_table, compare_strictly);
    CGAL_optimisation_assertion(
      M_new->number_of_columns() == M_new->number_of_rows());
  
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
  int* t_new = new int[M_new->number_of_rows() + 1];
  t_new[M_new->number_of_rows()] = M_new->number_of_columns();
  
  if ( M_new->number_of_rows() == 1)
    // recursion anchor:
    // we have just one element ==> no choice
    t_new[0] = 0;
  else
    monotone_matrix_search( *M_new, t_new);
  

  // and conquer
  // -----------
  
  int j( 0);       // actual index in t
  int j_new( 0);   // actual index in t_new
  do {
    // even row ==> we know
    *(t+j) = reduction_table[t_new[j_new++]];
    if ( ++j >= M.number_of_rows())
      break;
  
    // odd row
    // search *(t+j) between *(t+j-1) and t_new[j_new]:
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

} // monotone_matrix_search( M, t)


template < class Matrix,
           class RandomAccessIterator,
           class Compare_strictly >
void
_reduce_matrix(
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
} // _reduce_matrix( M, t)
} //namespace CGAL

#endif // ! (CGAL_MONOTONE_MATRIX_SEARCH_H)
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------
