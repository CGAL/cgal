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
// file          : include/CGAL/optimal_convex_partiton_2_test.C
// package       : $CGAL_Package: Partition_2 1.0 (27 Jul 2000) $
// chapter       : Planar Polygon Partitioning
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
// maintainer    : Susan Hert <hert@mpi-sb.mpg.de>
//
// coordinator   : MPI (Susan Hert <hert@mpi-sb.mpg.de>)
//
// implementation: Testing of convex partitioning functions
// ============================================================================

#include "kernel_include.h"
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_traits_2.h>
#include <CGAL/partition_2.h>
#include <list>
#include <cassert>

#include "kernel_def.h"

typedef K::Point_2                         Point_2;
typedef CGAL::Polygon_traits_2<K>          Traits;
typedef std::list<Point_2>                 Container;
typedef CGAL::Polygon_2<Traits, Container> Polygon_2;

#include "convex_test_polys.h"

void test_optimal_convex()
{
   Polygon_2              polygon;
   std::list<Polygon_2>   partition_polys;

   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_monotone_convex(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));

   assert (partition_polys.size() == 1 && 
           partition_polys.front().size() == polygon.size());
   assert(CGAL::is_convex_2(partition_polys.front().vertices_begin(), 
                            partition_polys.front().vertices_end()));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_convex_w_collinear_points(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex_w_collinear_points(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                   polygon.vertices_end(),
                                   std::back_inserter(partition_polys));

   partition_polys.clear();
   polygon.erase(polygon.vertices_begin(), polygon.vertices_end());
   make_nonconvex(polygon);
   CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));
}

int main(void)
{
   test_optimal_convex();

   return 0;
}
