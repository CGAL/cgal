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
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
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
  CGAL::convex_hull_2( V.begin(), V.end(), 
                       std::back_inserter(CH),
                       CGAL::Convex_hull_leda_traits_2() ); 
  return 0;
}
#endif // CGAL_USE_LEDA
