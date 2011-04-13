// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 03
//
// file          : ch_timing_2.C
// package       : Convex_hull (3.3)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : convex_hull_2.lw
// revision      : 3.3
// revision_date : 03 Aug 2000
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_TIMING_2_C
#define CGAL_CH_TIMING_2_C

#ifndef CGAL_CH_TIMING_2_H
#include <CGAL/ch_timing_2.h>
#endif // CGAL_CH_TIMING_2_H

CGAL_BEGIN_NAMESPACE
template <class ForwardIterator, class OutputIterator, class Traits>
void
ch_timing( ForwardIterator first, ForwardIterator last,
           OutputIterator result,
           int iterations, 
           const Traits& ch_traits)
{
#ifndef CGAL_CFG_NO_NAMESPACE
  using std::cout;
  using std::endl;
#endif // CGAL_CFG_NO_NAMESPACE

  int i;
  CGAL::Timer  Clck;
  double       delta_t;

  cout << endl;
  OutputIterator  restart = result;

  Clck.start();

  for (i=0; i < iterations; i++)
  {
      result = restart;
      ch_akl_toussaint( first, last , result, ch_traits);
  }
  Clck.stop();
  delta_t = Clck.time();
  Clck.reset();

  cout << "ch_akl_toussaint:         " << delta_t << endl;

  Clck.start();

  for (i=0; i < iterations; i++)
  {
      result = restart;
      ch_eddy( first, last , result, ch_traits);
  }
  Clck.stop();
  delta_t = Clck.time();
  Clck.reset();

  cout << "ch_eddy:                  " << delta_t << endl;

  Clck.start();

  for (i=0; i < iterations; i++)
  {
      result = restart;
      ch_bykat( first, last , result, ch_traits);
  }
  Clck.stop();
  delta_t = Clck.time();
  Clck.reset();

  cout << "ch_bykat                  " << delta_t << endl;

  Clck.start();

  for (i=0; i < iterations; i++)
  {
      result = restart;
      ch_bykat_with_threshold( first, last , result, ch_traits);
  }
  Clck.stop();
  delta_t = Clck.time();
  Clck.reset();

  cout << "ch_bykat_with_threshold:  " << delta_t << endl;

  Clck.start();

  for (i=0; i < iterations; i++)
  {
      result = restart;
      ch_graham_andrew( first, last , result, ch_traits);
  }
  Clck.stop();
  delta_t = Clck.time();
  Clck.reset();

  cout << "ch_graham_andrew:         " << delta_t << endl;

  Clck.start();

  for (i=0; i < iterations; i++)
  {
      result = restart;
      ch_jarvis( first, last , result, ch_traits);
  }
  Clck.stop();
  delta_t = Clck.time();
  Clck.reset();

  cout << "ch_jarvis:                " << delta_t << endl;

}

CGAL_END_NAMESPACE

#endif // CGAL_CH_TIMING_2_C

