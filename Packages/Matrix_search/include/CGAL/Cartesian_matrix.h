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
// file          : Cartesian_matrix.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// A Representation for Cartesian Matrices
// ============================================================================

#if ! (CGAL_CARTESIAN_MATRIX_H)
#define CGAL_CARTESIAN_MATRIX_H 1

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#include <CGAL/optimisation_assertions.h>
#endif // CGAL_OPTIMISATION_ASSERTIONS_H

template < class Operation,
           class RandomAccessIC_row,
           class RandomAccessIC_column >
class CGAL_Cartesian_matrix {
public:
  typedef typename Operation::result_type           Value;
  typedef typename Operation::first_argument_type   RowValue;
  typedef typename Operation::second_argument_type  ColumnValue;

  /*
  CGAL_Cartesian_matrix( Operation o = Operation())
  : op( o)
  {}
  */

  CGAL_Cartesian_matrix( RandomAccessIC_row r_f,
                         RandomAccessIC_row r_l,
                         RandomAccessIC_column c_f,
                         RandomAccessIC_column c_l)
  : row_vec( r_f),
    column_vec( c_f),
    n_rows( r_l - r_f),
    n_columns( c_l - c_f)
  {}

  CGAL_Cartesian_matrix( RandomAccessIC_row r_f,
                         RandomAccessIC_row r_l,
                         RandomAccessIC_column c_f,
                         RandomAccessIC_column c_l,
                         const Operation& o)
  : row_vec( r_f),
    column_vec( c_f),
    n_rows( r_l - r_f),
    n_columns( c_l - c_f),
    op( o)
  {}

  int
  number_of_rows() const
  { return n_rows; }

  int
  number_of_columns() const
  { return n_columns; }

  Value
  operator()( int r, int c) const
  {
    CGAL_optimisation_precondition( r >= 0 && r < number_of_rows());
    CGAL_optimisation_precondition( c >= 0 && c < number_of_columns());
    return op( row_vec[r], column_vec[c]);
  }

protected:
  RandomAccessIC_row     row_vec;
  RandomAccessIC_column  column_vec;
  int                    n_rows;
  int                    n_columns;
  Operation              op;
}; // class CGAL_Cartesian_matrix< ... >

template < class Operation,
           class RandomAccessIC_row,
           class RandomAccessIC_column >
inline
CGAL_Cartesian_matrix< Operation,
                       RandomAccessIC_row,
                       RandomAccessIC_column >
CGAL_cartesian_matrix( RandomAccessIC_row r_f,
                       RandomAccessIC_row r_l,
                       RandomAccessIC_column c_f,
                       RandomAccessIC_column c_l,
                       const Operation& o)
{
  return
  CGAL_Cartesian_matrix< Operation,
                         RandomAccessIC_row,
                         RandomAccessIC_column >
  ( r_f, r_l, c_f, c_l, o);
}

#endif // ! (CGAL_CARTESIAN_MATRIX_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

