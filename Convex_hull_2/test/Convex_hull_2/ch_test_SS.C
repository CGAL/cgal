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
// file          : ch_test_SS.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
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
  CGAL::Simple_cartesian<leda_rational>                ch_S_rational;
  std::cout << "SimpleCartesian<rational>:     ";
  CGAL::ch__batch_test( ch_S_rational );
#else
  CGAL::Simple_cartesian<CGAL::Quotient<CGAL::Gmpz> >  ch_S_Qgmp;
  std::cout << "SimpleCartesian<Quotient<Gmpz> > >:     ";
  CGAL::ch__batch_test( ch_S_Qgmp );
#endif // CGAL_USE_LEDA
  CGAL::Simple_cartesian<double>                       ch_S_double;
  std::cout << "SimpleCartesian<double>:     ";
  CGAL::ch__batch_test( ch_S_double );
  return 0;
}
