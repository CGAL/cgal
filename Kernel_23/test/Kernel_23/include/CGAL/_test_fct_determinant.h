// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Stefan Schirra


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

  // determinant()

  assert( CGAL::determinant<NT>( zero, one, \
                                   one,  zero) == - one );

  assert( CGAL::determinant<NT>( one,  one, \
                                   zero, -one) == - one );

  assert( CGAL::determinant<NT>( one,  zero, one, \
                                   zero, one,  one, \
                                  -one,  -one, zero ) == \
          one * CGAL::determinant<NT>( one, one, \
                                        -one, zero ) \
        - one * CGAL::determinant<NT>( zero, one, \
                                         one, one ) );

  assert( CGAL::determinant<NT>( one,  zero, one, -one, \
                                   zero, one,  one, -one, \
                                   one,  two,  zero, one, \
                                  -one,  -one, zero, zero ) == \
          one * CGAL::determinant<NT>( zero, one, -one, \
                                         one, one, -one, \
                                         two, zero, one ) \
        - one * CGAL::determinant<NT>( one, one, -one, \
                                         zero, one, -one, \
                                         one, zero, one ) );

  assert( CGAL::determinant<NT>( one,  zero, one, -one, zero, \
                                   zero, one,  one, -one, zero, \
                                   one,  two,  zero, one, two, \
                                   one, one,  zero, -one, zero, \
                                  -one,  -one, zero, zero, one ) == \
          two * CGAL::determinant<NT>( one,  zero, one,  -one, \
                                         zero, one,  one,  -one, \
                                         one,  one,  zero, -one, \
                                        -one, -one,  zero, zero) \
        + one * CGAL::determinant<NT>( one,  zero, one,  -one, \
                                         zero, one,  one,  -one, \
                                         one,  two,  zero,  one, \
                                         one,  one,  zero, -one) );

  assert( CGAL::determinant<NT>( one,  zero, zero, zero, zero, zero, \
                                   zero, one,  zero, one,  -one, zero, \
                                   zero, zero, one,  one, -one, zero, \
                                   zero, one,  two,  zero, one, two, \
                                   zero, one, one,  zero, -one, zero, \
                                   zero, -one,  -one, zero, zero, one ) == \
          CGAL::determinant<NT>( one,  zero, one, -one, zero, \
                                   zero, one,  one, -one, zero, \
                                   one,  two,  zero, one, two, \
                                   one, one,  zero, -one, zero, \
                                  -one,  -one, zero, zero, one ) );

  // sign_of_determinantDxD

  assert( CGAL::sign_of_determinant<NT>( zero, one, \
                                        one,  zero) == CGAL::NEGATIVE );

  assert( CGAL::sign_of_determinant<NT>( one,  one, \
                                        zero, -one) == CGAL::NEGATIVE );

  assert( CGAL::sign_of_determinant<NT>( one,  one, \
                                        zero, one) == CGAL::POSITIVE );

  assert( CGAL::sign_of_determinant<NT>( zero, -one, \
                                        zero, one) == CGAL::ZERO );


  assert( CGAL::sign_of_determinant<NT>( one,  zero, one, -one, \
                                        zero, one,  one, -one, \
                                        zero, two,  zero, one, \
                                        zero, -one, zero, zero ) == \
          CGAL::sign_of_determinant<NT>( one,  one, -one, \
                                        two,  zero, one, \
                                        -one, zero, zero ) );

  assert( CGAL::sign_of_determinant<NT>( one, zero, zero, zero, zero, zero, \
                                        zero, one,  zero, one,  -one, zero, \
                                        zero, zero, one,  one, -one, zero, \
                                        zero, one,  two,  zero, one, two, \
                                        zero, one, one,  zero, -one, zero, \
                                        zero, -one,  -one, zero, zero, one )==\
          CGAL::sign_of_determinant<NT>( one,  zero, one, -one, zero, \
                                        zero, one,  one, -one, zero, \
                                        one,  two,  zero, one, two, \
                                        one, one,  zero, -one, zero, \
                                       -one,  -one, zero, zero, one ) );


  return true;
}
#endif // CGAL__TEST_DETERMINANT_C
