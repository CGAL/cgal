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
// file          : ch__test.h
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL_CH__TEST_H
#define CGAL_CH__TEST_H

#include <CGAL/convex_hull_2.h>
#include <CGAL/convexity_check_2.h>
#include <CGAL/ch_akl_toussaint.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/ch_eddy.h>
#include <CGAL/ch_bykat.h>
#include <CGAL/ch_jarvis.h>


namespace CGAL {

enum ch_Algorithm{ ch_JARVIS,
                        ch_GRAHAM_ANDREW,
                        ch_EDDY,
                        ch_BYKAT,
                        ch_BYKAT_WITH_THRESHOLD,
                        ch_LOWER_UPPER,
                        ch_AKL_TOUSSAINT,
                        ch_ALL,
                        ch_DEFAULT };

enum ch_Check_status{ ch_CHECK_ALL,
                           ch_CHECK_CONVEXITY,
                           ch_CHECK_CONTAINEMENT,
                           ch_NO_CHECK };

template <class InputIterator, class Traits>
bool
ch__test(InputIterator first, InputIterator last,
              const Traits& ch_traits,
              ch_Algorithm alg,
              ch_Check_status check_level);

template <class InputIterator, class Traits>
bool
ch__test(InputIterator first, InputIterator last,
              const Traits& ch_traits,
              ch_Algorithm alg);

template <class InputIterator, class Traits>
bool
ch__test(InputIterator first, InputIterator last,
              const Traits& ch_traits);


} //namespace CGAL

#include <CGAL/ch__test_impl.h>

#endif // CGAL_CH__TEST_H
