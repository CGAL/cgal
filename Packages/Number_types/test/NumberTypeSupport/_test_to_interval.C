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
// file          : _test_to_interval.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <CGAL/Quotient.h>
#include <CGAL/_test_to_interval.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_real.h>
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif

#ifdef CGAL_USE_CLN
#include <CGAL/CLN/cl_integer.h>
#include <CGAL/CLN/cl_rational.h>
#endif

int
main()
{
  bool ok =
     test_to_interval( double())
  && test_to_interval( float())
  && test_to_interval( short())
  && test_to_interval( int())
  && test_to_interval( long())
  // && test_to_interval( long long())
  && test_to_interval( CGAL::Quotient< double>() )
  && test_to_interval( CGAL::Quotient< float>() )

#ifdef CGAL_USE_LEDA
  && test_to_interval( leda_bigfloat() )
  && test_to_interval( leda_integer() )
  && test_to_interval( leda_rational() )
  && test_to_interval( leda_real() )
  && test_to_interval( CGAL::Quotient< leda_bigfloat>() )
  && test_to_interval( CGAL::Quotient< leda_integer>() )
  && test_to_interval( CGAL::Quotient< leda_real>() )
#endif // CGAL_USE_LEDA

#ifdef CGAL_USE_GMP
  && test_to_interval( CGAL::Gmpz() )
  && test_to_interval( CGAL::Quotient< CGAL::Gmpz>() )
#endif // CGAL_USE_GMP

#ifdef CGAL_USE_CLN
  && test_to_interval( cl_I() )
  && test_to_interval( cl_RA() )
  && test_to_interval( CGAL::Quotient< cl_I>() )
#endif // CGAL_USE_CLN
  ;

  return ok ? 0 : 1;
}
