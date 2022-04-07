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
// file          : ch_test_CH.C
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ============================================================================

#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/convex_hull_traits_2.h>
#include <CGAL/convex_hull_constructive_traits_2.h>

#include <fstream>

#include <list>
#include <vector>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#endif// CGAL_USE_LEDA

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#endif

#include <CGAL/_test_fct_ch_I_2.h>

int
main()
{
#ifdef CGAL_USE_LEDA
  CGAL::Convex_hull_constructive_traits_2< CGAL::Homogeneous<leda_integer> >
                                                                 cch_H_integer;
  std::cout << "Homogeneous<integer>: C   ";
  CGAL::ch__batch_test( cch_H_integer );
#endif

#ifdef CGAL_USE_GMP
  CGAL::Convex_hull_constructive_traits_2< CGAL::Homogeneous<CGAL::Gmpz> >
                                                                 cch_H_gmp;
  std::cout << "Homogeneous<gmp>: C   ";
  CGAL::ch__batch_test( cch_H_gmp );
#endif

  CGAL::Convex_hull_constructive_traits_2< CGAL::Homogeneous<double> >
                                                                 cch_H_double;
  std::cout << "Homogeneous<double>: C ";
  CGAL::ch__batch_test( cch_H_double );
  return 0;
}
