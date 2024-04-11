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
#include <CGAL/ch_melkman.h>

#include <CGAL/Polygon_2_algorithms.h>

#include <array>
#include <iterator>
#include <set>
#include <string>
#include <vector>

namespace CGAL {

enum ch_Algorithm
{
  ch_AKL_TOUSSAINT = 0,
  ch_BYKAT,
  ch_BYKAT_WITH_THRESHOLD,
  ch_EDDY,
  ch_JARVIS,
  ch_GRAHAM_ANDREW,
  ch_LOWER_UPPER,
  ch_MELKMAN,
  ch_ALL,
  ch_DEFAULT
};

static std::array<std::string, 10> algorithm_names = {{ "AKL TOUSSAINT",
                                                        "BYKAT",
                                                        "BYKAT WITH THRESHOLD",
                                                        "EDDY",
                                                        "JARVIS",
                                                        "GRAHAM ANDREW",
                                                        "LOWER UPPER",
                                                        "MELKMAN",
                                                        "ALL",
                                                        "DEFAULT"
                                                      }};

enum ch_Check_status
{
  ch_CHECK_ALL,
  ch_CHECK_CONVEXITY,
  ch_CHECK_CONTAINEMENT,
  ch_NO_CHECK
};

// Also compare the results
template <typename InputIterator, typename Traits>
bool test_all_and_compare(InputIterator first, InputIterator beyond,
                          const Traits& ch_traits,
                          ch_Check_status check_level = ch_CHECK_ALL)
{
  using Point_2 = typename Traits::Point_2;

  bool call_melkman = is_simple_2(first, beyond, ch_traits);

  std::cout << "Input:";
  InputIterator first_o = first;
  for(; first_o!=beyond;++first_o)
    std::cout << " " << first_o->x() << " " << first_o->y() << " 0";
  std::cout << std::endl;

  if(!ch__test(first, beyond, ch_traits, ch_AKL_TOUSSAINT, check_level))
    return false;
  if(!ch__test(first, beyond, ch_traits, ch_BYKAT, check_level))
    return false;
  if(!ch__test(first, beyond, ch_traits, ch_BYKAT_WITH_THRESHOLD, check_level))
    return false;
  if(!ch__test(first, beyond, ch_traits, ch_EDDY, check_level))
    return false;
  if(!ch__test(first, beyond, ch_traits, ch_JARVIS, check_level))
    return false;
  if(!ch__test(first, beyond, ch_traits, ch_GRAHAM_ANDREW, check_level))
    return false;
  if(!ch__test(first, beyond, ch_traits, ch_LOWER_UPPER, check_level))
    return false;
  if(call_melkman)
  {
    if(!ch__test(first, beyond, ch_traits, ch_MELKMAN, check_level))
      return false;
  }

  std::set<Point_2> V0;
  ch_akl_toussaint(first, beyond, std::inserter(V0, V0.end()), ch_traits);

  {
    std::set<Point_2> V0_other;
    ch_bykat(first, beyond, std::inserter(V0_other, V0_other.end()), ch_traits);
    assert(V0 == V0_other);
  }

  {
    std::set<Point_2> V0_other;
    ch_bykat_with_threshold(first, beyond, std::inserter(V0_other, V0_other.end()), ch_traits);
    assert(V0 == V0_other);
  }

  {
    std::set<Point_2> V0_other;
    ch_eddy(first, beyond, std::inserter(V0_other, V0_other.end()), ch_traits);
    assert(V0 == V0_other);
  }

  {
    std::set<Point_2> V0_other;
    ch_jarvis(first, beyond, std::inserter(V0_other, V0_other.end()), ch_traits);
    assert(V0 == V0_other);
  }

  {
    std::set<Point_2> V0_other;
    ch_graham_andrew(first, beyond, std::inserter(V0_other, V0_other.end()), ch_traits);
    assert(V0 == V0_other);
  }

  {
    std::set<Point_2> V0_other;
    lower_hull_points_2(first, beyond, std::inserter(V0_other, V0_other.end()), ch_traits);
    upper_hull_points_2(first, beyond, std::inserter(V0_other, V0_other.end()), ch_traits);
    assert(V0 == V0_other);
  }

  if(call_melkman)
  {
    std::cout << "Simple polyline, proceed with melkman check" << std::endl;
    std::set<Point_2> V0_other;
    ch_melkman(first, beyond, std::inserter(V0_other, V0_other.end()), ch_traits);

    assert(V0 == V0_other);
  }
  else
  {
    std::cout << "Cannot call melkman on non-simple polyline" << std::endl;
  }

  return true;
}

template <class InputIterator, class Traits>
bool
ch__test(InputIterator first, InputIterator beyond,
         const Traits& ch_traits,
         ch_Algorithm alg = ch_ALL,
         ch_Check_status check_level = ch_CHECK_ALL)
{
  using Point_2 = typename Traits::Point_2;
  std::vector<Point_2> VI(first, beyond);
  std::vector<Point_2> VO;

  std::cout << "Algorithm: " << algorithm_names[std::size_t(alg)] << std::endl;
  switch(alg)
  {
    case ch_AKL_TOUSSAINT:
      ch_akl_toussaint(first, beyond, std::back_inserter(VO), ch_traits);
      break;
    case ch_BYKAT:
      ch_bykat(first, beyond, std::back_inserter(VO), ch_traits);
      break;
    case ch_BYKAT_WITH_THRESHOLD:
      ch_bykat_with_threshold(first, beyond, std::back_inserter(VO), ch_traits);
      break;
    case ch_EDDY:
      ch_eddy(first, beyond, std::back_inserter(VO), ch_traits);
      break;
    case ch_JARVIS:
      ch_jarvis(first, beyond, std::back_inserter(VO), ch_traits);
      break;
    case ch_GRAHAM_ANDREW:
      ch_graham_andrew(first, beyond, std::back_inserter(VO), ch_traits);
      break;
    case ch_LOWER_UPPER:
      lower_hull_points_2(first, beyond, std::back_inserter(VO), ch_traits);
      upper_hull_points_2(first, beyond, std::back_inserter(VO), ch_traits);
      break;
    case ch_MELKMAN:
      ch_melkman(first, beyond, std::back_inserter(VO), ch_traits);
      break;
    case ch_ALL:
      return test_all_and_compare(first, beyond, ch_traits, check_level);
    case ch_DEFAULT:
    default:
      convex_hull_2(first, beyond, std::back_inserter(VO), ch_traits);
      break;
  }

  switch(check_level)
  {
    case ch_CHECK_CONVEXITY:
      return is_ccw_strongly_convex_2(VO.begin(), VO.end(), ch_traits);
    case ch_CHECK_CONTAINEMENT:
      return ch_brute_force_check_2(first, beyond, VO.begin(), VO.end(), ch_traits);
    case ch_NO_CHECK:
      return true;
    case ch_CHECK_ALL:
    default:
      return is_ccw_strongly_convex_2(VO.begin(), VO.end(), ch_traits) &&
             ch_brute_force_check_2(first, beyond, VO.begin(), VO.end(), ch_traits);
  }
}

} // namespace CGAL

#endif // CGAL_CH__TEST_H
