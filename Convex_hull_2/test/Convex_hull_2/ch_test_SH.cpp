// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// 
// 
// 
//
// ----------------------------------------------------------------------------
// release       :
// release_date  :
//
// file          : ch_test_SH.C
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ============================================================================

#include <CGAL/Homogeneous.h>
#include <CGAL/convex_hull_traits_2.h>
#include <CGAL/convex_hull_constructive_traits_2.h>

#include <fstream>

#include <list>
#include <vector>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif
#include <CGAL/_test_fct_ch_I_2.h>

int
main()
{
#ifdef CGAL_USE_LEDA
  CGAL::Homogeneous<leda_integer>    ch_H_integer;
  std::cout << "Homogeneous<integer>:     ";
  CGAL::ch__batch_test( ch_H_integer );
#endif

#ifdef CGAL_USE_GMP
  CGAL::Homogeneous<CGAL::Gmpz>      ch_H_gmp;
  std::cout << "Homogeneous<gmp>:     ";
  CGAL::ch__batch_test( ch_H_gmp );
#endif

  CGAL::Homogeneous<double>           ch_H_double;
  std::cout << "Homogeneous<double>:   ";
  CGAL::ch__batch_test( ch_H_double );
  return 0;
}
