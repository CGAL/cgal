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
// file          : rectangular_p_center_2_example.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2-4-Centering Axis-Parallel 2D-Rectangles - example program
// ============================================================================

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_POINT_GENERATORS_2_H
#include <CGAL/point_generators_2.h>
#endif // CGAL_POINT_GENERATORS_2_H
#ifndef CGAL_RECTANGULAR_P_CENTER_2_H
#include <CGAL/rectangular_p_center_2.h>
#endif // CGAL_RECTANGULAR_P_CENTER_2_H
#ifndef CGAL_COPY_N_H
#include <CGAL/copy_n.h>
#endif // CGAL_COPY_N_H
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

typedef double                            FT;
typedef CGAL_Cartesian< FT >              R;
typedef CGAL_Point_2< R >                 Point_2;
typedef vector< Point_2 >                 Cont;
typedef CGAL_Random_points_in_square_2<
  Point_2,
  CGAL_Creator_uniform_2< FT, Point_2 > >
Point_generator;
typedef CGAL_Ostream_iterator< Point_2, ostream >
  Ostream_iterator_point;

int main() {

  int number_of_points( 10);
  int p( 2);
  Ostream_iterator_point cout_ip( cout);
  CGAL_set_pretty_mode( cout);

  Cont points;
  CGAL_copy_n( Point_generator( 1),
               number_of_points,
               back_inserter( points));
  cout << "Generated Point Set:" << endl;
  copy( points.begin(), points.end(), cout_ip);

  Cont centers;
  FT p_radius;
  CGAL_rectangular_p_center_2(
    points.begin(),
    points.end(),
    back_inserter( centers),
    p_radius,
    p);
  cout << "\n\n" << p << "-centers:" << endl;
  copy( centers.begin(), centers.end(), cout_ip);
  cout << "\n\n" << p << "-radius = " << p_radius << endl;

} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

