// ============================================================================
//
// Copyright (c) 2000 The GALIA Consortium
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
// file          : test/Generator/random_poly_test.C
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Random Simple Polygons: Test Program
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Polygon_2.h>
#include <list>

typedef CGAL::Cartesian< double >                           CR;
typedef CGAL::Point_2< CR >                                 CPoint_2;
typedef std::list<CPoint_2>                                 CContainer;
typedef CGAL::Polygon_2<CR, CContainer>                     CPolygon_2;
typedef CGAL::Creator_uniform_2<double, CPoint_2>           CCreator;
typedef CGAL::Random_points_in_square_2<CPoint_2, CCreator> CPoint_generator;

typedef CGAL::Homogeneous< double >                          HR;
typedef CGAL::Point_2< HR >                                  HPoint_2;
typedef std::list<HPoint_2>                                  HContainer;
typedef CGAL::Polygon_2<HR, HContainer>                      HPolygon_2;
typedef CGAL::Creator_uniform_2<double, HPoint_2>            HCreator;
typedef CGAL::Random_points_in_square_2<HPoint_2, HCreator>  HPoint_generator;

int main() {

  CPolygon_2 polygon1;
  int n = 50;

  // create a polygon 
  CGAL::random_polygon_2(n, std::back_inserter(polygon1), 
                         CPoint_generator(0.5));

  // make sure it is simple
  if (! polygon1.is_simple())
  {
     std::cerr << "ERROR: polygon is not simple." << std::endl;
     return 1;
  }

  HPolygon_2 polygon2;

  // create a polygon 
  CGAL::random_polygon_2(n, std::back_inserter(polygon2), 
                         HPoint_generator(0.5));

  // make sure it is simple
  if (! polygon1.is_simple())
  {
     std::cerr << "ERROR: polygon is not simple." << std::endl;
     return 1;
  }
  return 0;
}

// EOF
