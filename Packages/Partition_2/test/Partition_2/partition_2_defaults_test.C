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
// file          : include/CGAL/partition_2_default_test.C
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
// implementation: testing of use of default traits classes for partitioning
// ============================================================================

#include "kernel_include.h"
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <list>
#include <cassert>

#include "kernel_def.h"

typedef CGAL::Partition_traits_2<K>                       Traits;
typedef Traits::Point_2                                   Point_2;
typedef Traits::Polygon_2                                 Polygon_2;
typedef std::list<Polygon_2>                              Polygon_list;


void make_polygon(Polygon_2& polygon)
{
   polygon.push_back(Point_2(227,423));
   polygon.push_back(Point_2(123,364));
   polygon.push_back(Point_2(129,254));
   polygon.push_back(Point_2(230,285));
   polygon.push_back(Point_2(231,128));
   polygon.push_back(Point_2(387,205));
   polygon.push_back(Point_2(417,331));
   polygon.push_back(Point_2(319,225));
   polygon.push_back(Point_2(268,293));
   polygon.push_back(Point_2(367,399));
   polygon.push_back(Point_2(298,418));
   polygon.push_back(Point_2(196,326));
}


int main( int argc, char** argv )
{
   Polygon_2    polygon;
   Polygon_list partition_polys;

   make_polygon(polygon);

   // this calls y_monotone_partition_is_valid_2 for postcondition checking
   // with a traits class specified, which calls partition_is_valid_2
   // with Is_y_monotone and supplied traits
   CGAL::y_monotone_partition_2(polygon.vertices_begin(), 
                                polygon.vertices_end(),
                                std::back_inserter(partition_polys));


   // check y-monotone and nonoverlapping with default traits
   CGAL::y_monotone_partition_is_valid_2(polygon.vertices_begin(),
                                         polygon.vertices_end(),
                                         partition_polys.begin(),
                                         partition_polys.end());
   // check y-monotone with traits supplied
   assert(CGAL::is_y_monotone_2((*partition_polys.begin()).vertices_begin(),
                                (*partition_polys.begin()).vertices_end(),
                                Traits()));
   // check y-monotone with default traits 
   assert(CGAL::is_y_monotone_2((*partition_polys.begin()).vertices_begin(),
                                (*partition_polys.begin()).vertices_end()));


   // checks for overlapping polygons using default traits
   assert(CGAL::partition_is_valid_2(polygon.vertices_begin(),
                                     polygon.vertices_end(),
                                     partition_polys.begin(),
                                     partition_polys.end()));

   typedef Traits::Is_y_monotone_2      Is_y_monotone_2;
   CGAL::Partition_is_valid_traits_2<Traits, Is_y_monotone_2>  validity_traits;

   // checks for overlapping polygons using supplied traits
   assert(CGAL::partition_is_valid_2(polygon.vertices_begin(),
                                     polygon.vertices_end(),
                                     partition_polys.begin(),
                                     partition_polys.end(),
                                     validity_traits));

   partition_polys.clear();
   // this calls convex_partition_is_valid_2 for postcondition checking
   // with a traits class specified, which calls partition_is_valid_2
   // with Is_convex_2 and supplied traits
   CGAL::approx_convex_partition_2(polygon.vertices_begin(), 
                                   polygon.vertices_end(),
                                   std::back_inserter(partition_polys));

   // check convex and nonoverlapping with default traits
   CGAL::convex_partition_is_valid_2(polygon.vertices_begin(),
                                     polygon.vertices_end(),
                                     partition_polys.begin(),
                                     partition_polys.end());

   return 0;
}
