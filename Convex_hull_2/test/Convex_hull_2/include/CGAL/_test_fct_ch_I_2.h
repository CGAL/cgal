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
// file          : _test_fct_ch_I_2.h
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ============================================================================

#ifndef _TEST_FCT_CH_I_2_H
#define _TEST_FCT_CH_I_2_H

#include <CGAL/ch__test.h>

#include <CGAL/Random.h>

#include <cassert>

namespace CGAL {

template <class Traits>
void ch__batch_test_simple(const Traits& chI)
{
  using Point_2 = typename Traits::Point_2;

  std::cout << std::endl << " == Test Simple ==" << std::endl;

  std::vector<Point_2> points;
  points.push_back(Point_2( -10,   0, 1));
  points.push_back(Point_2( -15,  -1, 1));
  points.push_back(Point_2(  -7, -10, 1));
  points.push_back(Point_2(  -9, -20, 1));
  points.push_back(Point_2(  -3, -20, 1));
  points.push_back(Point_2(  74,  0, 1));
  points.push_back(Point_2(   1,  1, 1));
  points.push_back(Point_2(   0, 100, 1));

  assert(ch__test(points.begin(), points.end(), chI));
}

template <class Traits>
void ch__batch_test_cocircular(const Traits& chI)
{
  using Point_2 = typename Traits::Point_2;

  std::cout << std::endl << " == Test Cocircular ==" << std::endl;

  std::vector<Point_2> cocircular_points;
  cocircular_points.push_back(Point_2(  39,  80,  89));
  cocircular_points.push_back(Point_2( 180, 299, 349));
  cocircular_points.push_back(Point_2(  -3,  -4,   5));
  cocircular_points.push_back(Point_2(-651, 260, 701));
  cocircular_points.push_back(Point_2( 180, -19, 181));
  cocircular_points.push_back(Point_2(-153, 104, 185));
  cocircular_points.push_back(Point_2(-247, -96, 265));
  cocircular_points.push_back(Point_2( -32, 255, 257));
  cocircular_points.push_back(Point_2(  45, -28,  53));
  cocircular_points.push_back(Point_2( -12, -35,  37));

  assert(!ch_brute_force_check_2(cocircular_points.begin(), cocircular_points.end(),
                                 cocircular_points.begin(), cocircular_points.end(), chI));
  assert(ch_brute_force_check_2(cocircular_points.begin(), cocircular_points.begin(),
                                cocircular_points.begin(), cocircular_points.end(), chI));
  assert(ch__test(cocircular_points.begin(), cocircular_points.end(), chI, ch_ALL, ch_CHECK_CONVEXITY));

  std::vector<Point_2> extreme_points;
  convex_hull_2(cocircular_points.begin(), cocircular_points.end(),
                std::back_inserter(extreme_points), chI);

  assert(is_ccw_strongly_convex_2(extreme_points.begin(), extreme_points.begin(), chI));
  assert(is_cw_strongly_convex_2(extreme_points.rend(), extreme_points.rend(), chI));
  assert(is_ccw_strongly_convex_2(extreme_points.begin(), extreme_points.begin() + 1, chI));
  assert(is_cw_strongly_convex_2(extreme_points.begin(), extreme_points.begin() + 1, chI));
  assert(is_ccw_strongly_convex_2(extreme_points.begin(), extreme_points.end(), chI));
  assert(is_cw_strongly_convex_2(extreme_points.rbegin(), extreme_points.rend(), chI));
}

template <class Traits>
void ch__batch_test_cocircular_ordered(const Traits& chI)
{
  using Point_2 = typename Traits::Point_2;

  std::cout << std::endl << " == Test Cocircular (ordered) ==" << std::endl;

  std::vector<Point_2> cocircular_points;
  cocircular_points.push_back(Point_2( 180, -19, 181));
  cocircular_points.push_back(Point_2(  45, -28,  53));
  cocircular_points.push_back(Point_2( -12, -35,  37));
  cocircular_points.push_back(Point_2(  -3,  -4,   5));
  cocircular_points.push_back(Point_2(-247, -96, 265));
  cocircular_points.push_back(Point_2(-651, 260, 701));
  cocircular_points.push_back(Point_2(-153, 104, 185));
  cocircular_points.push_back(Point_2( -32, 255, 257));
  cocircular_points.push_back(Point_2(  39,  80,  89));
  cocircular_points.push_back(Point_2( 180, 299, 349));

  assert(!ch_brute_force_check_2(cocircular_points.begin(), cocircular_points.end(),
                                 cocircular_points.begin(), cocircular_points.end(), chI));
  assert(ch_brute_force_check_2(cocircular_points.begin(), cocircular_points.begin(),
                                cocircular_points.begin(), cocircular_points.end(), chI));
  assert(ch__test(cocircular_points.begin(), cocircular_points.end(), chI, ch_ALL, ch_CHECK_CONVEXITY));

  std::vector<Point_2> extreme_points;
  convex_hull_2(cocircular_points.begin(), cocircular_points.end(),
                std::back_inserter(extreme_points), chI);

  assert(is_ccw_strongly_convex_2(extreme_points.begin(), extreme_points.begin(), chI));
  assert(is_cw_strongly_convex_2(extreme_points.rend(), extreme_points.rend(), chI));
  assert(is_ccw_strongly_convex_2(extreme_points.begin(), extreme_points.begin() + 1, chI));
  assert(is_cw_strongly_convex_2(extreme_points.begin(), extreme_points.begin() + 1, chI));
  assert(is_ccw_strongly_convex_2(extreme_points.begin(), extreme_points.end(), chI));
  assert(is_cw_strongly_convex_2(extreme_points.rbegin(), extreme_points.rend(), chI));
}

template <class Traits>
void ch__batch_test_collinear(const Traits& chI)
{
  using Point_2 = typename Traits::Point_2;

  std::cout << std::endl << " == Test Collinear ==" << std::endl;

  std::vector<Point_2> collinear_points;
  collinear_points.push_back(Point_2( 16,  20, 1));
  collinear_points.push_back(Point_2( 46,  40, 1));
  collinear_points.push_back(Point_2( 76,  60, 1));
  collinear_points.push_back(Point_2(106,  80, 1));
  collinear_points.push_back(Point_2(-14,   0, 1));
  collinear_points.push_back(Point_2(136, 100, 1));

  assert(ch__test(collinear_points.begin(), collinear_points.end(), chI));
}

template <class Traits>
void ch__batch_test_multiple(const Traits& chI)
{
  using Point_2 = typename Traits::Point_2;

  std::cout << std::endl << " == Test Multiple ==" << std::endl;

  std::vector<Point_2> multiple_points;
  multiple_points.push_back(Point_2(17, 80, 1));
  multiple_points.push_back(Point_2(17, 80, 1));
  multiple_points.push_back(Point_2(17, 80, 1));
  multiple_points.push_back(Point_2(17, 80, 1));

  assert(ch_brute_force_check_2(multiple_points.begin(), multiple_points.end(),
                                multiple_points.begin(), multiple_points.begin() + 1, chI));
  assert(is_ccw_strongly_convex_2(multiple_points.begin(), multiple_points.begin(), chI));

  assert(ch__test(multiple_points.begin(), multiple_points.end(), chI));
  assert(ch__test(multiple_points.begin() + 2, multiple_points.begin() + 3, chI));
}

template <class Traits>
void ch__batch_test_iso_rectangle(const Traits& chI)
{
  using Point_2 = typename Traits::Point_2;

  std::cout << std::endl << " == Test Iso Rectangle ==" << std::endl;

  std::vector<Point_2> iso_rectangle_points;
  iso_rectangle_points.push_back(Point_2( 15,   0,  1));
  iso_rectangle_points.push_back(Point_2( 45,   0,  1));
  iso_rectangle_points.push_back(Point_2( 70,   0, 10));
  iso_rectangle_points.push_back(Point_2( 12,   0,  1));
  iso_rectangle_points.push_back(Point_2( 56, 118,  1));
  iso_rectangle_points.push_back(Point_2( 27, 118,  1));
  iso_rectangle_points.push_back(Point_2( 56, 118,  1));
  iso_rectangle_points.push_back(Point_2(112, 118,  1));
  iso_rectangle_points.push_back(Point_2(  0,   9,  1));
  iso_rectangle_points.push_back(Point_2(  0,  78,  1));
  iso_rectangle_points.push_back(Point_2(  0,  16,  1));
  iso_rectangle_points.push_back(Point_2(  0,  77,  1));
  iso_rectangle_points.push_back(Point_2(150,  56,  1));
  iso_rectangle_points.push_back(Point_2(150,  57,  1));
  iso_rectangle_points.push_back(Point_2(150,  58,  1));
  iso_rectangle_points.push_back(Point_2(150,  58,  1));

  assert(ch__test(iso_rectangle_points.begin(), iso_rectangle_points.end(), chI));
  assert(ch__test(iso_rectangle_points.begin(), iso_rectangle_points.begin()+3, chI));
  assert(ch__test(iso_rectangle_points.begin()+4, iso_rectangle_points.begin()+7, chI));
  assert(ch__test(iso_rectangle_points.begin()+5, iso_rectangle_points.begin()+5, chI));
}

// https://github.com/CGAL/cgal/issues/6723
template <class Traits>
void ch__batch_test_issue_6723(const Traits& chI)
{
  using Point_2 = typename Traits::Point_2;

  std::cout << std::endl << " == Test issue_6723 ==" << std::endl;

  std::vector<Point_2> points = { Point_2( 4, 2, 1),
                                  Point_2( 0, 0, 1),
                                  Point_2(10, 0, 1),
                                  Point_2( 3, 1, 1),
                                  Point_2( 3,-1, 1),
                                  Point_2( 2, 2, 1),
                                  Point_2( 5, 2, 1),
                                  Point_2( 9, 2, 1),
                                  Point_2( 3, 2, 1) };

  assert(ch__test(points.begin(), points.end(), chI));
}

template <class Traits>
void ch__batch_test_random(const Traits& chI)
{
  using Point_2 = typename Traits::Point_2;

  std::cout << std::endl << " == Test Random ==" << std::endl;

  CGAL::Random rnd;
  std::cout << "Random seed = " << rnd.get_seed() << std::endl;

  std::vector<Point_2> points;

  const int h = rnd.get_int(1, 100);
  for(int i=0; i<10; ++i)
  {
    const int x = rnd.get_int(-100, 100);
    const int y = rnd.get_int(-100, 100);
    points.emplace_back(x, y, h);
  }

  assert(ch__test(points.begin(), points.end(), chI));
}

template <class Traits>
bool ch__batch_test(const Traits& chI)
{
  std::cout << "Testing Convex Hulls" << std::endl;;

  ch__batch_test_simple(chI);
  ch__batch_test_cocircular(chI);
  ch__batch_test_cocircular_ordered(chI); // so melkman gets called
  ch__batch_test_collinear(chI);
  ch__batch_test_multiple(chI);
  ch__batch_test_iso_rectangle(chI);
  ch__batch_test_issue_6723(chI);
  ch__batch_test_random(chI);

  std::cout << "done" << std::endl;
  return true;
}

} // namespace CGAL

#endif // _TEST_FCT_CH_I_2_H
