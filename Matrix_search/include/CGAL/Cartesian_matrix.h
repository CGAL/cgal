// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_CARTESIAN_MATRIX_H
#define CGAL_CARTESIAN_MATRIX_H 1

#include <CGAL/license/Matrix_search.h>


#include <CGAL/basic.h>
#include <CGAL/Optimisation/assertions.h>

namespace CGAL {

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
    n_rows(static_cast<int>(r_l - r_f)),
    n_columns(static_cast<int>(c_l - c_f))
  {}

  Cartesian_matrix(RandomAccessIC_row r_f,
                   RandomAccessIC_row r_l,
                   RandomAccessIC_column c_f,
                   RandomAccessIC_column c_l,
                   const Operation& o)
  : row_vec(r_f),
    column_vec(c_f),
    n_rows(static_cast<int>(r_l - r_f)),
    n_columns(static_cast<int>(c_l - c_f)),
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

} //namespace CGAL

#endif // ! (CGAL_CARTESIAN_MATRIX_H)
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------
