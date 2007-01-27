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
// file          : include/CGAL/partition_2_default_test_C.C
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
// implementation: testing of use of default traits classes for partitioning
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <list>
#include <cassert>

typedef CGAL::Cartesian<double>       K;
typedef CGAL::Partition_traits_2<K>   Traits;
typedef Traits::Point_2               Point_2;
typedef Traits::Polygon_2             Polygon_2;
typedef std::list<Polygon_2>          Polygon_list;

#include "test_defaults.h"

int main( )
{
   test_defaults();

   return 0;
}
