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
// file          : Sorted_matrix_search_traits_adaptor.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : fjsearch.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Frederickson-Johnson matrix search: traits class adaptor
// ============================================================================

#if ! (CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H)
#define CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H 1

CGAL_BEGIN_NAMESPACE

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

CGAL_END_NAMESPACE

#endif // ! (CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

