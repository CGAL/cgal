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

#ifndef CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H
#define CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H 1

#include <CGAL/license/Matrix_search.h>


#include <functional>

namespace CGAL {

template < class FeasibilityTest_, class Matrix_ >
class Sorted_matrix_search_traits_adaptor {
public:
  typedef FeasibilityTest_         FeasibilityTest;
  typedef Matrix_                  Matrix;
  typedef typename Matrix::Value   Value;
  typedef std::less< Value >       Compare_strictly;
  typedef std::less_equal< Value > Compare_non_strictly;

  Sorted_matrix_search_traits_adaptor(FeasibilityTest ft)
  : ft_( ft)
  {}

  Compare_strictly
  compare_strictly() const
  { return Compare_strictly(); }

  Compare_non_strictly
  compare_non_strictly() const
  { return Compare_non_strictly(); }

  bool
  is_feasible(Value a)
  { return ft_( a); }

protected:
  FeasibilityTest ft_;
};

//!!! with iterator traits we replace const Matrix&
// by an iterator with value type Matrix
template < class FeasibilityTest, class Matrix >
Sorted_matrix_search_traits_adaptor<
  FeasibilityTest, Matrix >
sorted_matrix_search_traits_adaptor(FeasibilityTest f, const Matrix&)
{
  typedef Sorted_matrix_search_traits_adaptor<
    FeasibilityTest, Matrix > Traits;
  return Traits(f);
} // sorted_matrix_search_traits_adaptor( ... )

} //namespace CGAL

#endif // ! (CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H)
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------
