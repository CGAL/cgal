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


#include <CGAL/config.h>
#include <CGAL/float.h>
#include <CGAL/number_utils.h>
#include <cassert>
#include <cmath>
#include <CGAL/IEEE_754_unions.h>

#if defined(BOOST_MSVC)
#  pragma warning(disable:4723)
#endif

bool
_test_valid_finite( const float& )
{
  float d0 = 0.0;
  float d1 = 1.0;
  float d2;

  assert( CGAL::is_valid( d0) );
  assert( CGAL::is_valid( d1) );

  assert( CGAL_NTS is_finite( d0) );
  assert( CGAL_NTS is_finite( d1) );

  if ( CGAL::is_valid( d1/d0 - d1/d0 ))
  { d2 = d1/d0 - d1/d0; show( reinterpret_cast<IEEE_754_float*>(&d2)); }
  if ( CGAL_NTS is_finite( d1/d0 ))
  { d2 = d1/d0; show( reinterpret_cast<IEEE_754_float*>(&d2)); }

  assert( CGAL::is_valid( d1/d0 ));
  assert( !CGAL::is_valid( d1/d0 - d1/d0 ));
  assert( !CGAL_NTS is_finite( d1/d0 ));

  return true;
}

int
main()
{
  _test_valid_finite( (float)1.0  );

  return 0;
}
