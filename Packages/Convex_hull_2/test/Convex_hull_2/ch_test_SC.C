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
// file          : ch_test_SC.C
// source        : convex_hull_2.lw
// revision      : 3.3
// revision_date : 03 Aug 2000
// author(s)     : Stefan Schirra
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================

#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>

#include <CGAL/convex_hull_cartesian_double_traits_2.h>

#include <CGAL/convex_hull_traits_2.h>

#include <CGAL/convex_hull_constructive_traits_2.h>

#include <fstream>

#include <deque>
#include <list>
#include <vector>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#else
#include <CGAL/Gmpz.h>
#endif // CGAL_USE_LEDA
#include <CGAL/_test_fct_ch_I_2.h>

int
main()
{
#ifdef CGAL_USE_LEDA
  CGAL::Cartesian<leda_rational>                ch_C_rational;
  std::cout << "Cartesian<rational>:     ";
  CGAL::ch__batch_test( ch_C_rational );
#else
  CGAL::Cartesian<CGAL::Quotient<CGAL::Gmpz> >  ch_C_Qgmp;
  std::cout << "Cartesian<Quotient<Gmpz> > >:     ";
  CGAL::ch__batch_test( ch_C_Qgmp );
#endif // CGAL_USE_LEDA
  CGAL::Cartesian<double>                       ch_C_double;
  std::cout << "Cartesian<double>:     ";
  CGAL::ch__batch_test( ch_C_double );
  return 0;
}
