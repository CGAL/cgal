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
// file          : ch_test_L.C
// source        : convex_hull_2.lw
// revision      : 3.3
// revision_date : 03 Aug 2000
// author(s)     : Stefan Schirra
// maintainer    : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL_USE_LEDA
int main() { return 0; }
#else
#include <CGAL/basic.h>
#include <vector>
#include <LEDA/point.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_leda_traits_2.h>

typedef leda_point  Point2;

int
main() 
{
  std::vector< Point2 > V;
  std::ifstream F("data/CD500");
  CGAL::set_ascii_mode( F );
  std::istream_iterator< Point2>  in_start( F );
  std::istream_iterator< Point2>  in_end;
  std::copy( in_start, in_end , std::back_inserter(V) );
  std::vector< Point2 > CH;
  CGAL::convex_hull_points_2( V.begin(), V.end(), 
                              std::back_inserter(CH),
                              CGAL::convex_hull_leda_traits_2() ); 
  return 0;
}
#endif // CGAL_USE_LEDA
