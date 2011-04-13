// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        :
// file          : _test_fct_determinant.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_DETERMINANT_C
#define CGAL__TEST_DETERMINANT_C

#include <CGAL/determinant.h>
#include <CGAL/predicates/sign_of_determinant.h>

template <class NT>
bool
_test_fct_determinant(const NT&)
{
  NT zero= NT(0);
  NT one = NT(1);
  NT two = NT(2);

  // detDxD_by_formula

  assert( CGAL::det2x2_by_formula( zero, one, \
                                   one,  zero) == - one );

  assert( CGAL::det2x2_by_formula( one,  one, \
                                   zero, -one) == - one );

  assert( CGAL::det3x3_by_formula( one,  zero, one, \
                                   zero, one,  one, \
                                  -one,  -one, zero ) == \
          one * CGAL::det2x2_by_formula( one, one, \
                                        -one, zero ) \
        - one * CGAL::det2x2_by_formula( zero, one, \
                                         one, one ) );

  assert( CGAL::det4x4_by_formula( one,  zero, one, -one, \
                                   zero, one,  one, -one, \
                                   one,  two,  zero, one, \
                                  -one,  -one, zero, zero ) == \
          one * CGAL::det3x3_by_formula( zero, one, -one, \
                                         one, one, -one, \
                                         two, zero, one ) \
        - one * CGAL::det3x3_by_formula( one, one, -one, \
                                         zero, one, -one, \
                                         one, zero, one ) );

  assert( CGAL::det5x5_by_formula( one,  zero, one, -one, zero, \
                                   zero, one,  one, -one, zero, \
                                   one,  two,  zero, one, two, \
                                   one, one,  zero, -one, zero, \
                                  -one,  -one, zero, zero, one ) == \
          two * CGAL::det4x4_by_formula( one,  zero, one,  -one, \
                                         zero, one,  one,  -one, \
                                         one,  one,  zero, -one, \
                                        -one, -one,  zero, zero) \
        + one * CGAL::det4x4_by_formula( one,  zero, one,  -one, \
                                         zero, one,  one,  -one, \
                                         one,  two,  zero,  one, \
                                         one,  one,  zero, -one) );

  assert( CGAL::det6x6_by_formula( one,  zero, zero, zero, zero, zero, \
                                   zero, one,  zero, one,  -one, zero, \
                                   zero, zero, one,  one, -one, zero, \
                                   zero, one,  two,  zero, one, two, \
                                   zero, one, one,  zero, -one, zero, \
                                   zero, -one,  -one, zero, zero, one ) == \
          CGAL::det5x5_by_formula( one,  zero, one, -one, zero, \
                                   zero, one,  one, -one, zero, \
                                   one,  two,  zero, one, two, \
                                   one, one,  zero, -one, zero, \
                                  -one,  -one, zero, zero, one ) );

  // sign_of_determinantDxD

  assert( CGAL::sign_of_determinant2x2( zero, one, \
                                        one,  zero) == CGAL::NEGATIVE );

  assert( CGAL::sign_of_determinant2x2( one,  one, \
                                        zero, -one) == CGAL::NEGATIVE );

  assert( CGAL::sign_of_determinant2x2( one,  one, \
                                        zero, one) == CGAL::POSITIVE );

  assert( CGAL::sign_of_determinant2x2( zero, -one, \
                                        zero, one) == CGAL::ZERO );


  assert( CGAL::sign_of_determinant4x4( one,  zero, one, -one, \
                                        zero, one,  one, -one, \
                                        zero, two,  zero, one, \
                                        zero, -one, zero, zero ) == \
          CGAL::sign_of_determinant3x3( one,  one, -one, \
                                        two,  zero, one, \
                                        -one, zero, zero ) );

  assert( CGAL::sign_of_determinant6x6( one,  zero, zero, zero, zero, zero, \
                                        zero, one,  zero, one,  -one, zero, \
                                        zero, zero, one,  one, -one, zero, \
                                        zero, one,  two,  zero, one, two, \
                                        zero, one, one,  zero, -one, zero, \
                                        zero, -one,  -one, zero, zero, one )==\
          CGAL::sign_of_determinant5x5( one,  zero, one, -one, zero, \
                                        zero, one,  one, -one, zero, \
                                        one,  two,  zero, one, two, \
                                        one, one,  zero, -one, zero, \
                                       -one,  -one, zero, zero, one ) );


  return true;
}
#endif // CGAL__TEST_DETERMINANT_C
