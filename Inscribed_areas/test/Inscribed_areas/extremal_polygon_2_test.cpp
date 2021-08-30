// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/extremal_polygon_2.h>
#include <functional>
#include <vector>

typedef double                                    FT;

typedef CGAL::Simple_cartesian<FT>                Kernel;

typedef Kernel::Point_2                           Point;
typedef std::vector<Point>                        Cont;
typedef CGAL::Polygon_2<Kernel>                   Polygon_2;
typedef CGAL::Random_points_in_square_2<Point>    Generator;

template < class RandomAccessIC,
           class OutputIterator >
OutputIterator
brute_force_area_3( RandomAccessIC b,
                    RandomAccessIC e,
                    OutputIterator o)
{
  RandomAccessIC i1( b);
  RandomAccessIC m1, m2, m3;
  FT a_m( -1);
  do {
    RandomAccessIC i2( i1);
    do {
      RandomAccessIC i3( i2);
      do {
        FT a(
          CGAL_NTS abs(
            (*i1).x() * ( (*i2).y() - (*i3).y()) +
            (*i2).x() * ( (*i3).y() - (*i1).y()) +
            (*i3).x() * ( (*i1).y() - (*i2).y())));
        if ( a > a_m) {
          a_m = a;
          m1 = i1;
          m2 = i2;
          m3 = i3;
        }
      } while ( ++i3 != e);
    } while ( ++i2 != e);
  } while ( ++i1 != e);

  *o++ = *m1;
  *o++ = *m2;
  *o++ = *m3;
  return o;

} // brute_force_3( b, e, o)
template < class RandomAccessIC,
           class OutputIterator >
OutputIterator
brute_force_area_4( RandomAccessIC b,
                    RandomAccessIC e,
                    OutputIterator o)
{
  RandomAccessIC i1( b);
  RandomAccessIC m1, m2, m3, m4;
  FT a_m( -1);
  do {
    RandomAccessIC i2( i1);
    do {
      RandomAccessIC i3( i2);
      do {
        RandomAccessIC i4( i3);
        do {
          FT a(
            CGAL_NTS abs(
              (*i1).x() * ( (*i4).y() - (*i3).y()) +
              (*i4).x() * ( (*i3).y() - (*i1).y()) +
              (*i3).x() * ( (*i1).y() - (*i4).y())) +
            CGAL_NTS abs(
              (*i1).x() * ( (*i2).y() - (*i3).y()) +
              (*i2).x() * ( (*i3).y() - (*i1).y()) +
              (*i3).x() * ( (*i1).y() - (*i2).y())));
          if ( a > a_m) {
            a_m = a;
            m1 = i1;
            m2 = i2;
            m3 = i3;
            m4 = i4;
          }
        } while ( ++i4 != e);
      } while ( ++i3 != e);
    } while ( ++i2 != e);
  } while ( ++i1 != e);

  *o++ = *m1;
  *o++ = *m2;
  *o++ = *m3;
  *o++ = *m4;
  return o;

} // brute_force_4( b, e, o)


int main() {

  int number_of_points[] = { 20, 51, 102, 500 };
  int k[] = { 3, 7, 12, 27 };
  int j;

  for (int n = 0; n < 4; ++n) {
    Polygon_2 p;
    CGAL::random_convex_set_2(number_of_points[n],
                              std::back_inserter(p),
                              Generator(1));
    assert(p.is_convex());

    for ( j = 0; j < 4; ++j) {
      // maximum area:
      Cont k_gon;
      k_gon.reserve( k[j]);

      CGAL::maximum_area_inscribed_k_gon_2(
        p.vertices_begin(),  p.vertices_end(),
        k[j],                std::back_inserter(k_gon));

    } // for ( j = 0; j < 4; ++j)

    for ( j = 0; j < 4; ++j) {
      // maximum perimeter:
      Cont k_gon;
      k_gon.reserve( k[j]);
      CGAL::maximum_perimeter_inscribed_k_gon_2(
        p.vertices_begin(),  p.vertices_end(),
        k[j],                std::back_inserter( k_gon));
    } // for ( j = 0; j < 4; ++j)

  } // for ( int n( 0); n < 4; ++n)

  return 0;
} // int main()
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

