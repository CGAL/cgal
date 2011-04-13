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
// file          : Sorted_matrix_search_traits_adaptor.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : fjsearch.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Frederickson-Johnson matrix search: traits class adaptor
// ============================================================================

#if ! (CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H)
#define CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H 1

template < class _FeasibilityTest, class _Matrix >
class CGAL_Sorted_matrix_search_traits_adaptor {
public:
  typedef _FeasibilityTest         FeasibilityTest;
  typedef _Matrix                  Matrix;
  typedef typename _Matrix::Value  Value;
  typedef less< Value >            Compare_strictly;
  typedef less_equal< Value >      Compare_non_strictly;

  CGAL_Sorted_matrix_search_traits_adaptor(
    const FeasibilityTest& ft)
  : _ft( ft)
  {}

  Compare_strictly
  compare_strictly() const
  { return Compare_strictly(); }

  Compare_non_strictly
  compare_non_strictly() const
  { return Compare_non_strictly(); }

  bool
  is_feasible( Value a)
  { return _ft( a); }

protected:
  FeasibilityTest _ft;
};

//!!! with iterator traits we replace const Matrix&
// by an iterator with value type Matrix
template < class FeasibilityTest, class Matrix >
CGAL_Sorted_matrix_search_traits_adaptor<
  FeasibilityTest, Matrix >
CGAL_sorted_matrix_search_traits_adaptor(
  const FeasibilityTest& f, const Matrix&)
{
  typedef CGAL_Sorted_matrix_search_traits_adaptor<
    FeasibilityTest, Matrix > Traits;
  return Traits( f);
} // CGAL_sorted_matrix_search_traits_adaptor( ... )

#endif // ! (CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

