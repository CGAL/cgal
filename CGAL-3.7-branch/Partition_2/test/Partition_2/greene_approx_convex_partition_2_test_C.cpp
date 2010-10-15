// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/greene_approx_convex_partiton_2_test_C.C
// package       : $CGAL_Package: Partition_2 1.0 (27 Jul 2000) $
// chapter       : Planar Polygon Partitioning
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Testing of convex partitioning functions
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/partition_2.h>
#include <list>
#include <cassert>

typedef CGAL::Cartesian<double>       K;
typedef K::Point_2                    Point_2;
typedef std::list<Point_2>            Container;
typedef CGAL::Polygon_2<K, Container> Polygon_2;

#include "convex_test_polys.h"

#include "test_greene_approx_convex.h"

int main(void)
{
   test_greene_approx_convex();

   return 0;
}
