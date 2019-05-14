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
// file          : test/Generator/rcs_test.C
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// source        : src/rcs/rcs.aw
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Random Convex Point Sets: Test Program
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <vector>
#include <iostream>
#include <cstdlib>

using std::vector;
using std::back_inserter;
using CGAL::Cartesian;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::set_pretty_mode;
using CGAL::random_convex_set_2;


int
main( )
{
  typedef Cartesian< double >                            R;
  typedef CGAL::Point_2< R >                             Point_2;
  typedef vector< Point_2 >                              Cont;
  typedef CGAL::Polygon_2< R, Cont >                     Polygon_2;
  typedef Creator_uniform_2< double, Point_2 >           Creator;
  typedef Random_points_in_square_2< Point_2, Creator >  Point_generator;
  
  // this is not initialized on MIPSPRO:
  set_pretty_mode( std::cout);
  set_pretty_mode( std::cerr);

  // define polygon:
  Polygon_2 p;
  int n( 1000);

  // build random n-gon:
  random_convex_set_2( n, back_inserter( p), Point_generator( 1));

  // check convexity:
  if ( ! p.is_convex()) {
    std::cerr << "ERROR: polygon is not convex." << std::endl;
    return 1;
  }

  return 0;
} // int main( )

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

