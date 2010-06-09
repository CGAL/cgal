// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of an example program for CGAL.  This example
// program may be used, distributed and modified without limitation.

#ifndef CGAL_CH_TIMING_2_H
#define CGAL_CH_TIMING_2_H

#include <CGAL/Timer.h>

namespace CGAL {

template <class ForwardIterator, class OutputIterator, class Traits>
void
ch_timing( ForwardIterator first, ForwardIterator last,
           OutputIterator result,
           int iterations, 
           const Traits& ch_traits);

} //namespace CGAL

#include <CGAL/ch_timing_2_impl.h>

#endif // CGAL_CH_TIMING_2_H
