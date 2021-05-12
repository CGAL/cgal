// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of an example program for CGAL.  This example
// program may be used, distributed and modified without limitation.

#ifndef CGAL_CH_TIMING_2_IMPL_H
#define CGAL_CH_TIMING_2_IMPL_H

#include <CGAL/ch_timing_2.h>

namespace CGAL {

template <class ForwardIterator, class OutputIterator, class Traits>
void
ch_timing( ForwardIterator first, ForwardIterator last,
           OutputIterator result,
           int iterations,
           const Traits& ch_traits)
{
  using std::cout;
  using std::endl;

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

} //namespace CGAL

#endif // CGAL_CH_TIMING_2_IMPL_H
