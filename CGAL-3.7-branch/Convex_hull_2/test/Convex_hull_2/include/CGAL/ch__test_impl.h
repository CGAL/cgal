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
// file          : ch__test.C
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL_CH__TEST_IMPL_H
#define CGAL_CH__TEST_IMPL_H

#include <CGAL/ch__test.h>

namespace CGAL {
template <class InputIterator, class Traits>
bool
ch__test(InputIterator first, InputIterator last, const Traits& ch_traits)
{
  ch_Algorithm    alg         = ch_ALL;
  ch_Check_status check_level = ch_CHECK_ALL;
  typedef typename Traits::Point_2   Point_2;
  std::vector< Point_2 >     VI (first, last);
  std::vector< Point_2 >     VO;
  typedef typename std::vector< Point_2 >::iterator  V_iter;
  V_iter VIfirst = VI.begin();
  V_iter VIlast  = VI.end();
  switch (alg) 
  {
      case ch_JARVIS:  
                ch_jarvis(VIfirst, VIlast, std::back_inserter(VO), ch_traits); 
                break;
      case ch_GRAHAM_ANDREW:
                ch_graham_andrew(VIfirst, VIlast, std::back_inserter(VO), 
                                      ch_traits);
                break;
      case ch_EDDY:
                ch_eddy(VIfirst, VIlast, std::back_inserter(VO), 
                             ch_traits);
                break;
      case ch_AKL_TOUSSAINT:
                ch_akl_toussaint( VIfirst, VIlast, std::back_inserter(VO),
                                       ch_traits);
                break;
      case ch_BYKAT:
                ch_bykat(VIfirst, VIlast, std::back_inserter(VO), 
                              ch_traits);
                break;
      case ch_BYKAT_WITH_THRESHOLD:
                ch_bykat_with_threshold(VIfirst, VIlast, std::back_inserter(VO), 
                                             ch_traits);
                break;
      case ch_LOWER_UPPER:
                lower_hull_points_2(VIfirst, VIlast, std::back_inserter(VO),
                                         ch_traits);
                upper_hull_points_2(VIfirst, VIlast, std::back_inserter(VO),
                                         ch_traits);
                break;
      case ch_ALL:
                return
                     ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_JARVIS, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_GRAHAM_ANDREW, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_EDDY, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_BYKAT, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_BYKAT_WITH_THRESHOLD, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_LOWER_UPPER, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_AKL_TOUSSAINT, check_level);
      case ch_DEFAULT:
      default:
                convex_hull_2( VIfirst, VIlast, std::back_inserter(VO), 
                               ch_traits);
                break;
  }

  switch (check_level)
  {
      case ch_CHECK_CONVEXITY:
                return is_ccw_strongly_convex_2( VO.begin(), VO.end(), 
                                                      ch_traits);
      case ch_CHECK_CONTAINEMENT:
                return ch_brute_force_check_2( VIfirst, VIlast,
                                                    VO.begin(), VO.end(), 
                                                    ch_traits);
      case ch_NO_CHECK:
                return true;
      case ch_CHECK_ALL:
      default:
                return is_ccw_strongly_convex_2( VO.begin(), VO.end(), 
                                                      ch_traits)
                    && ch_brute_force_check_2( VIfirst, VIlast,
                                                    VO.begin(), VO.end(), 
                                                    ch_traits);
  }
}

template <class InputIterator, class Traits>
bool
ch__test(InputIterator first, InputIterator last, const Traits& ch_traits, 
              ch_Algorithm alg )
{
  ch_Check_status check_level = ch_CHECK_ALL;
  typedef typename Traits::Point_2   Point_2;
  std::vector< Point_2 >     VI (first, last);
  std::vector< Point_2 >     VO;
  typedef typename std::vector< Point_2 >::iterator  V_iter;
  V_iter VIfirst = VI.begin();
  V_iter VIlast  = VI.end();
  switch (alg) 
  {
      case ch_JARVIS:  
                ch_jarvis(VIfirst, VIlast, std::back_inserter(VO), ch_traits); 
                break;
      case ch_GRAHAM_ANDREW:
                ch_graham_andrew(VIfirst, VIlast, std::back_inserter(VO), 
                                      ch_traits);
                break;
      case ch_EDDY:
                ch_eddy(VIfirst, VIlast, std::back_inserter(VO), 
                             ch_traits);
                break;
      case ch_AKL_TOUSSAINT:
                ch_akl_toussaint( VIfirst, VIlast, std::back_inserter(VO),
                                       ch_traits);
                break;
      case ch_BYKAT:
                ch_bykat(VIfirst, VIlast, std::back_inserter(VO), 
                              ch_traits);
                break;
      case ch_BYKAT_WITH_THRESHOLD:
                ch_bykat_with_threshold(VIfirst, VIlast, std::back_inserter(VO), 
                                             ch_traits);
                break;
      case ch_LOWER_UPPER:
                lower_hull_points_2(VIfirst, VIlast, std::back_inserter(VO),
                                         ch_traits);
                upper_hull_points_2(VIfirst, VIlast, std::back_inserter(VO),
                                         ch_traits);
                break;
      case ch_ALL:
                return
                     ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_JARVIS, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_GRAHAM_ANDREW, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_EDDY, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_BYKAT, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_BYKAT_WITH_THRESHOLD, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_LOWER_UPPER, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_AKL_TOUSSAINT, check_level);
      case ch_DEFAULT:
      default:
                ( VIfirst, VIlast, std::back_inserter(VO), 
                                           ch_traits);
                break;
  }

  switch (check_level)
  {
      case ch_CHECK_CONVEXITY:
                return is_ccw_strongly_convex_2( VO.begin(), VO.end(), 
                                                      ch_traits);
      case ch_CHECK_CONTAINEMENT:
                return ch_brute_force_check_2( VIfirst, VIlast,
                                                    VO.begin(), VO.end(), 
                                                    ch_traits);
      case ch_NO_CHECK:
                return true;
      case ch_CHECK_ALL:
      default:
                return is_ccw_strongly_convex_2( VO.begin(), VO.end(), 
                                                      ch_traits)
                    && ch_brute_force_check_2( VIfirst, VIlast,
                                                    VO.begin(), VO.end(), 
                                                    ch_traits);
  }
}

template <class InputIterator, class Traits>
bool
ch__test(InputIterator first, InputIterator last, const Traits& ch_traits, 
              ch_Algorithm alg, ch_Check_status check_level)
{
  typedef typename Traits::Point_2   Point_2;
  std::vector< Point_2 >     VI (first, last);
  std::vector< Point_2 >     VO;
  typedef typename std::vector< Point_2 >::iterator  V_iter;
  V_iter VIfirst = VI.begin();
  V_iter VIlast  = VI.end();
  switch (alg) 
  {
      case ch_JARVIS:  
                ch_jarvis(VIfirst, VIlast, std::back_inserter(VO), ch_traits); 
                break;
      case ch_GRAHAM_ANDREW:
                ch_graham_andrew(VIfirst, VIlast, std::back_inserter(VO), 
                                      ch_traits);
                break;
      case ch_EDDY:
                ch_eddy(VIfirst, VIlast, std::back_inserter(VO), 
                             ch_traits);
                break;
      case ch_AKL_TOUSSAINT:
                ch_akl_toussaint( VIfirst, VIlast, std::back_inserter(VO),
                                       ch_traits);
                break;
      case ch_BYKAT:
                ch_bykat(VIfirst, VIlast, std::back_inserter(VO), 
                              ch_traits);
                break;
      case ch_BYKAT_WITH_THRESHOLD:
                ch_bykat_with_threshold(VIfirst, VIlast, std::back_inserter(VO), 
                                             ch_traits);
                break;
      case ch_LOWER_UPPER:
                lower_hull_points_2(VIfirst, VIlast, std::back_inserter(VO),
                                         ch_traits);
                upper_hull_points_2(VIfirst, VIlast, std::back_inserter(VO),
                                         ch_traits);
                break;
      case ch_ALL:
                return
                     ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_JARVIS, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_GRAHAM_ANDREW, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_EDDY, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_BYKAT, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_BYKAT_WITH_THRESHOLD, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_LOWER_UPPER, check_level)
                  && ch__test(VIfirst, VIlast, ch_traits, 
                                   ch_AKL_TOUSSAINT, check_level);
      case ch_DEFAULT:
      default:
                convex_hull_2( VIfirst, VIlast, std::back_inserter(VO), 
                               ch_traits);
                break;
  }

  switch (check_level)
  {
      case ch_CHECK_CONVEXITY:
                return is_ccw_strongly_convex_2( VO.begin(), VO.end(), 
                                                      ch_traits);
      case ch_CHECK_CONTAINEMENT:
                return ch_brute_force_check_2( VIfirst, VIlast,
                                                    VO.begin(), VO.end(), 
                                                    ch_traits);
      case ch_NO_CHECK:
                return true;
      case ch_CHECK_ALL:
      default:
                return is_ccw_strongly_convex_2( VO.begin(), VO.end(), 
                                                      ch_traits)
                    && ch_brute_force_check_2( VIfirst, VIlast,
                                                    VO.begin(), VO.end(), 
                                                    ch_traits);
  }
}

} //namespace CGAL

#endif // CGAL_CH__TEST_IMPL_H
