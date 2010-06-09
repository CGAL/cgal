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
// file          : _test_fct_ch_I_2.C
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef _TEST_FCT_CH_I_2_IMPL_H
#define _TEST_FCT_CH_I_2_IMPL_H

#include <CGAL/_test_fct_ch_I_2.h>
#include <cassert>

namespace CGAL {

template <class Traits>
bool
ch__batch_test( const Traits& chI )
{
  typedef typename Traits::Point_2   Point_2;
 
  std::cout << "Testing ch";
  std::vector< Point_2 >  Cocircular_points;
  Cocircular_points.push_back( Point_2( 39, 80, 89 ));
  Cocircular_points.push_back( Point_2( 180, 299, 349 ));
  Cocircular_points.push_back( Point_2( -3, -4, 5 ));
  Cocircular_points.push_back( Point_2( -651, 260, 701 ));
  Cocircular_points.push_back( Point_2( 180, -19, 181 ));
  Cocircular_points.push_back( Point_2( -153, 104, 185 ));
  Cocircular_points.push_back( Point_2( -247, -96, 265 ));
  Cocircular_points.push_back( Point_2( -32, 255, 257 ));
  Cocircular_points.push_back( Point_2( 45, -28, 53 ));
  Cocircular_points.push_back( Point_2( -12, -35, 37 ));
  assert( ! ch_brute_force_check_2( \
              Cocircular_points.begin(), Cocircular_points.end(), \
              Cocircular_points.begin(), Cocircular_points.end(), chI ));
  assert( ch_brute_force_check_2( \
              Cocircular_points.begin(), Cocircular_points.begin(), \
              Cocircular_points.begin(), Cocircular_points.end(), chI ));
  assert( ch__test( Cocircular_points.begin(), \
                         Cocircular_points.end(), \
                         chI, ch_ALL, ch_CHECK_CONVEXITY ));
  std::vector< Point_2 > extreme_points;
  convex_hull_2(Cocircular_points.begin(), Cocircular_points.end(),
                std::back_inserter(extreme_points), chI );
  assert( is_ccw_strongly_convex_2( extreme_points.begin(), \
                                         extreme_points.begin(), \
                                         chI ));
  assert( is_cw_strongly_convex_2( extreme_points.rend(), 
                                        extreme_points.rend(), \
                                        chI ));
  assert( is_ccw_strongly_convex_2( extreme_points.begin(), \
                                extreme_points.begin() + 1,  chI ));
  assert( is_cw_strongly_convex_2( extreme_points.begin(), \
                                extreme_points.begin() + 1,  chI ));
  assert( is_ccw_strongly_convex_2( extreme_points.begin(), \
                                         extreme_points.end(), \
                                         chI ));
  assert( is_cw_strongly_convex_2( extreme_points.rbegin(), 
                                        extreme_points.rend(), \
                                        chI ));

  std::cout << '.';
  std::vector< Point_2 > Collinear_points;
  Collinear_points.push_back( Point_2( 16, 20, 1 ));
  Collinear_points.push_back( Point_2( 46, 40, 1 ));
  Collinear_points.push_back( Point_2( 76, 60, 1 ));
  Collinear_points.push_back( Point_2( 106, 80, 1 ));
  Collinear_points.push_back( Point_2( -14, 0 , 1));
  Collinear_points.push_back( Point_2( 136, 100, 1 ));
  assert( ch__test( Collinear_points.begin(), \
                         Collinear_points.end(), chI ));

  std::cout << '.';
  std::vector< Point_2 > Multiple_points;
  Multiple_points.push_back( Point_2( 17, 80, 1 ));
  Multiple_points.push_back( Point_2( 17, 80, 1 ));
  Multiple_points.push_back( Point_2( 17, 80, 1 ));
  Multiple_points.push_back( Point_2( 17, 80, 1 ));
  assert( ch_brute_force_check_2( \
               Multiple_points.begin(), Multiple_points.end(),\
               Multiple_points.begin(), Multiple_points.begin() + 1, chI ));
  assert( is_ccw_strongly_convex_2(Multiple_points.begin(), \
                                        Multiple_points.begin(), chI ));
  assert( ch__test( Multiple_points.begin(), \
                         Multiple_points.end(), chI ));
  assert( ch__test( Multiple_points.begin() + 2, \
                         Multiple_points.begin() + 3, chI ));

  std::cout << '.';
  std::vector< Point_2 > Iso_rectangle_points;
  Iso_rectangle_points.push_back( Point_2( 15, 0, 1 ));
  Iso_rectangle_points.push_back( Point_2( 45, 0, 1 ));
  Iso_rectangle_points.push_back( Point_2( 70, 0, 10 ));
  Iso_rectangle_points.push_back( Point_2( 12, 0, 1 ));
  Iso_rectangle_points.push_back( Point_2( 56, 118, 1 ));
  Iso_rectangle_points.push_back( Point_2( 27, 118, 1 ));
  Iso_rectangle_points.push_back( Point_2( 56, 118, 1 ));
  Iso_rectangle_points.push_back( Point_2( 112, 118, 1));
  Iso_rectangle_points.push_back( Point_2( 0, 9, 1));
  Iso_rectangle_points.push_back( Point_2( 0, 78, 1));
  Iso_rectangle_points.push_back( Point_2( 0, 16, 1));
  Iso_rectangle_points.push_back( Point_2( 0, 77, 1));
  Iso_rectangle_points.push_back( Point_2( 150, 56, 1));
  Iso_rectangle_points.push_back( Point_2( 150, 57, 1));
  Iso_rectangle_points.push_back( Point_2( 150, 58, 1));
  Iso_rectangle_points.push_back( Point_2( 150, 58, 1));
  assert( ch__test( Iso_rectangle_points.begin(), \
                         Iso_rectangle_points.end(), chI ));
  assert( ch__test( Iso_rectangle_points.begin(), \
                         Iso_rectangle_points.begin()+3, chI ));
  assert( ch__test( Iso_rectangle_points.begin()+4, \
                         Iso_rectangle_points.begin()+7, chI ));
  assert( ch__test( Iso_rectangle_points.begin()+5, \
                         Iso_rectangle_points.begin()+5, chI ));

  std::cout << "done" << std::endl;
  return true;
}

} //namespace CGAL

#endif // _TEST_FCT_CH_I_2_IMPL_H
