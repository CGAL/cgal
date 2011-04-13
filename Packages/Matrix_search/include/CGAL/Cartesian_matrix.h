// ============================================================================
//
// Copyright (c) 1998, 1999, 2000 The CGAL Consortium
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
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// A Representation for Cartesian Matrices
// ============================================================================

#if ! (CGAL_CARTESIAN_MATRIX_H)
#define CGAL_CARTESIAN_MATRIX_H 1

#include <CGAL/basic.h>
#include <CGAL/Optimisation/assertions.h>

CGAL_BEGIN_NAMESPACE

template < class Operation,
           class RandomAccessIC_row,
           class RandomAccessIC_column >
class Cartesian_matrix {
public:
  typedef typename Operation::result_type           Value;

  Cartesian_matrix(RandomAccessIC_row r_f,
                   RandomAccessIC_row r_l,
                   RandomAccessIC_column c_f,
                   RandomAccessIC_column c_l)
  : row_vec(r_f),
    column_vec(c_f),
    n_rows(r_l - r_f),
    n_columns(c_l - c_f)
  {}

  Cartesian_matrix(RandomAccessIC_row r_f,
                   RandomAccessIC_row r_l,
                   RandomAccessIC_column c_f,
                   RandomAccessIC_column c_l,
                   const Operation& o)
  : row_vec(r_f),
    column_vec(c_f),
    n_rows(r_l - r_f),
    n_columns(c_l - c_f),
    op(o)
  {}

  int
  number_of_rows() const
  { return n_rows; }

  int
  number_of_columns() const
  { return n_columns; }

  Value
  operator()(int r, int c) const
  {
    CGAL_optimisation_precondition(r >= 0 && r < number_of_rows());
    CGAL_optimisation_precondition(c >= 0 && c < number_of_columns());
    return op(row_vec[r], column_vec[c]);
  }

protected:
  RandomAccessIC_row     row_vec;
  RandomAccessIC_column  column_vec;
  int                    n_rows;
  int                    n_columns;
  Operation              op;
}; // class Cartesian_matrix< ... >

template < class Operation,
           class RandomAccessIC_row,
           class RandomAccessIC_column >
inline
Cartesian_matrix< Operation,
                       RandomAccessIC_row,
                       RandomAccessIC_column >
cartesian_matrix( RandomAccessIC_row r_f,
                       RandomAccessIC_row r_l,
                       RandomAccessIC_column c_f,
                       RandomAccessIC_column c_l,
                       const Operation& o)
{
  return
  Cartesian_matrix< Operation,
                         RandomAccessIC_row,
                         RandomAccessIC_column >
  ( r_f, r_l, c_f, c_l, o);
}

CGAL_END_NAMESPACE

#endif // ! (CGAL_CARTESIAN_MATRIX_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

