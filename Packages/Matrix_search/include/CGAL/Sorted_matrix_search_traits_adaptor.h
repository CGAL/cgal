#line 939 "fjsearch.aw"
#line 18 "code_formatting.awi"
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

#line 943 "fjsearch.aw"
#line 54 "code_formatting.awi"
#if ! (CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H)
#define CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H 1

#line 892 "fjsearch.aw"
#line 46 "code_formatting.awi"
CGAL_BEGIN_NAMESPACE
#line 893 "fjsearch.aw"

template < class _FeasibilityTest, class _Matrix >
class Sorted_matrix_search_traits_adaptor {
public:
  typedef _FeasibilityTest         FeasibilityTest;
  typedef _Matrix                  Matrix;
  typedef typename _Matrix::Value  Value;
  typedef std::less< Value >       Compare_strictly;
  typedef std::less_equal< Value > Compare_non_strictly;

  Sorted_matrix_search_traits_adaptor(FeasibilityTest ft)
  : _ft( ft)
  {}

  Compare_strictly
  compare_strictly() const
  { return Compare_strictly(); }

  Compare_non_strictly
  compare_non_strictly() const
  { return Compare_non_strictly(); }

  bool
  is_feasible(Value a)
  { return _ft( a); }

protected:
  FeasibilityTest _ft;
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

#line 50 "code_formatting.awi"
CGAL_END_NAMESPACE
#line 936 "fjsearch.aw"

#endif // ! (CGAL_SORTED_MATRIX_SEARCH_TRAITS_ADAPTOR_H)

#line 12 "code_formatting.awi"
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

