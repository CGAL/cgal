// ============================================================================
//
// Copyright (c) 1997-2004 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/QP_solver/NonstandardForm.C
// package       : $CGAL_Package: QP_engine $
// chapter       : Quadratic Programming Engine
//
// revision      : 3.0alpha
// revision_date : 2004/06
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: utility routines for QP's in nonstandard form
// ============================================================================

CGAL_BEGIN_NAMESPACE

// Looks in x_O_v_i which bound is present for variable i and returns
// the variable's value corresponding to this bound.
//
// Precondition: Is_in_standard_form is Tag_false.
//
// Notee: This routine is probably never called (but it is referenced
// sometimes in comments for explanation.)
template < class Rep_ >
typename QP_solver<Rep_>::ET
QP_solver<Rep_>::original_variable_value(int i) const
{
  CGAL_assertion(!check_tag(Is_in_standard_form()) && i<qp_n);
  switch (x_O_v_i[i]) {
  case UPPER:
    return qp_u[i];
    break;
  case ZERO:
    return et0;
    break;
  case LOWER:
  case FIXED:
    return qp_l[i];
    break;
  case BASIC:
    CGAL_qpe_assertion(false);
  }
  return et0; // dummy
}

// Returns the current value of a nonbasic original variable.
//
// Note: this routine is like original_variable_value() above with the
// two differences that (i) it also works if In_standard_form is
// Tag_true and (ii) that it makes the assertion check that i is not a
// basic variable.
//
// Precondition: x_O_v_i must be initialized as well as in_B. 
template < class Rep_ >
typename QP_solver<Rep_>::ET
QP_solver<Rep_>::nonbasic_original_variable_value(int i) const
{
  if (check_tag(Is_in_standard_form()))
    return et0;

  CGAL_assertion(!is_basic(i));
  return original_variable_value(i);
}

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
