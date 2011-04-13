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
// file          : test_kernel__.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/_test_misc.C>
#include <CGAL/_test_cls_quotient.C>
#include <CGAL/_test_fct_point_conversion.C>
#include <CGAL/_test_fct_determinant.C>
#include <CGAL/MP_Float.h>
#include <CGAL/Precise_numbers.h>
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#endif

int
main()
{
  std::cout << "Testing miscellanous" << std::endl;
  typedef   CGAL::Homogeneous<Precise_integer>                 H_Cls;
  typedef   CGAL::Cartesian<CGAL::Quotient<Precise_integer> >  C_Cls;
  typedef   CGAL::Simple_cartesian<CGAL::Quotient<Precise_integer> >  S_Cls;
  _test_cls_object( C_Cls() );
  _test_cls_object( H_Cls() );
  _test_cls_object( S_Cls() );
  _test_cls_quotient( int(1) );
  _test_cls_quotient( double(1.0) );
  _test_cls_quotient( long(1) );
  _test_cls_quotient( short(1) );
  _test_fct_determinant( Precise_integer(1) );
  _test_fct_point_conversion( int(1) );
  _test_fct_point_conversion( double(1.0) );
#ifdef CGAL_USE_GMP
  _test_cls_quotient( CGAL::Gmpz(1) );
  _test_fct_point_conversion( CGAL::Gmpz(1) );
#endif
#ifdef CGAL_USE_LEDA
  _test_cls_quotient( leda_integer(1) );
  _test_fct_point_conversion( leda_integer(1) );
#endif
  _test_cls_quotient( CGAL::MP_Float(1) );
  _test_fct_point_conversion( CGAL::MP_Float(1) );
  return 0;
}
