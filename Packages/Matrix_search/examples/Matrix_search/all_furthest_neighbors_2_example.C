// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : all_furthest_neighbors_2_example.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Example program: All Furthest Neighbors for a Convex Polygon
// ============================================================================

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_POLYGON_2_H
#include <CGAL/Polygon_2.h>
#endif // CGAL_POLYGON_2_H
#ifndef CGAL_POINT_GENERATORS_2_H
#include <CGAL/point_generators_2.h>
#endif // CGAL_POINT_GENERATORS_2_H
#ifndef CGAL_RANDOM_CONVEX_SET_2_H
#include <CGAL/random_convex_set_2.h>
#endif // CGAL_RANDOM_CONVEX_SET_2_H
#ifndef CGAL_ALL_FURTHEST_NEIGHBORS_2_H
#include <CGAL/all_furthest_neighbors_2.h>
#endif // CGAL_ALL_FURTHEST_NEIGHBORS_2_H
#ifndef CGAL_IO_OSTREAM_ITERATOR_H
#include <CGAL/IO/Ostream_iterator.h>
#endif // CGAL_IO_OSTREAM_ITERATOR_H
#ifndef CGAL_PROTECT_IOSTREAM_H
#include <iostream.h>
#define CGAL_PROTECT_IOSTREAM_H
#endif // CGAL_PROTECT_IOSTREAM_H
#ifndef CGAL_PROTECT_VECTOR_H
#include <vector.h>
#define CGAL_PROTECT_VECTOR_H
#endif // CGAL_PROTECT_VECTOR_H

typedef double                                 FT;
typedef CGAL_Cartesian< FT >                   R;
typedef CGAL_Point_2< R >                      Point_2;
typedef CGAL_Polygon_traits_2< R >             P_traits;
typedef vector< Point_2 >                      Point_cont;
typedef CGAL_Polygon_2< P_traits, Point_cont > Polygon_2;
typedef CGAL_Random_points_in_square_2<
  Point_2,
  CGAL_Creator_uniform_2< FT, Point_2 > >
Point_generator;
typedef CGAL_Ostream_iterator< int, ostream >  Ostream_iterator;

int
main()
{
  // generate random convex polygon:
  Polygon_2 p;
  CGAL_random_convex_set_2( 10,
                            back_inserter( p),
                            Point_generator( 1));

  // compute all furthest neighbors:
  CGAL_all_furthest_neighbors(
    p.vertices_begin(),
    p.vertices_end(),
    Ostream_iterator( cout));
  cout << endl;

} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

